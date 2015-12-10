//
//  arc-module.cpp
//  ndnrtc
//
//  Created by Takahiro Yoneda on 11/13/15.
//  Copyright 2013-2015 Regents of the University of California
//

#include <iostream>
#include <limits>
#include <cmath>
#include <cstring>

#include <boost/asio.hpp>
#include <boost/asio/steady_timer.hpp>
//#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
//#include <boost/ref.hpp>
#include <boost/system/error_code.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/tag.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

#include "arc-interface.h"
#include "arc-module.h"

using namespace ndnrtc;

const ArcModule::ArcIndicators ArcModule::ZeroIndicators = {0,0.,0.,0.,0.,0.,0.,0.,ConsumerPhaseInactive};

int ArcModule::initialize(IRateAdaptationModuleCallback* const callback,
                          const CodecMode& codecMode,
                          std::vector<ThreadEntry> mediaThreads)
{
    ArcTval tv;
    getNowTval(&tv);
    callback_ = callback;
    if (codecMode != CodecMode::CodecModeNormal) return -1;
    
    //registration rate control function to boost io service
    arcTimer_.expires_from_now(boost::chrono::milliseconds(ARC_INTERVAL));
    arcTimer_.async_wait(boost::bind(&ArcModule::autoRateControl, this));
    
    //initialize parameter
    numThread_ = 0;
    arcState_ = arcStateNormal;
    arcStateCheck_ = true;
    nextInterestPps_ = 0;
    countChallengePhase_ = 0;
    currThHist_ = nextThHist_ = NULL;
    arcCallTval_ = arcStartTval_ = arcCallThroughputTval_ = tv;
    sumDataSize_ = 0;

    for (auto itr = mediaThreads.begin(); itr != mediaThreads.end(); ++itr) {
        ThreadTable[numThread_] = *itr;
#ifdef ARC_DEBUG_INITIALIZE
        std::cout << "ArcModule initializing thread_id[" << ThreadTable[numThread_].id_ << "] bitrate[" <<  ThreadTable[numThread_].bitrate_ << "kbps] parity[" << ThreadTable[numThread_].parityRatio_ << "]" << std::endl;
#endif //ARC_DEBUG_INITIALIZE
        ++numThread_;
    }
    return 0;
}

void ArcModule::interestExpressed(const std::string &name,
                                  unsigned int threadId,
				  uint32_t interestNonce)
{
#ifdef ARC_DEBUG_SNDINTEREST
    std::cout << "ArcModule interestExpressed thread_id[" << threadId << "] name[" << name << "]" << std::endl;
#endif //ARC_DEBUG_SNDINTEREST
    if (consumerPhase_ == ConsumerPhaseAdjust)
        currThreadId_ = threadId;
    if (!(consumerPhase_ == ConsumerPhaseFetch || consumerPhase_ == ConsumerPhaseChallenge)) return;
    
    if (threadId == currThreadId_ && currThHist_ != NULL) {
        currThHist_->interestExpressed(name, threadId, interestNonce);
    } else if (threadId == nextThreadId_ && nextThHist_ != NULL) {
        nextThHist_->interestExpressed(name, threadId, interestNonce);
    }
    return;
}

void ArcModule::interestRetransmit(const std::string &name,
                                   unsigned int threadId)
{
#ifdef ARC_DEBUG_SNDINTEREST
    std::cout << "ArcModule interestRetransmit thread_id[" << threadId << "] name[" << name << "]" << std::endl;
#endif //ARC_DEBUG_SNDINTEREST
    if (!(consumerPhase_ == ConsumerPhaseFetch || consumerPhase_ == ConsumerPhaseChallenge)) return;
    
    if (threadId == currThreadId_ && currThHist_ != NULL) {
        currThHist_->interestRetransmit(name, threadId);
    } else if (threadId == nextThreadId_ && nextThHist_ != NULL) {
        nextThHist_->interestRetransmit(name, threadId);
    }
    return;
}

void ArcModule::dataReceived(const std::string &interestName,
                             const std::string &dataName,
                             unsigned int threadId,
                             unsigned int ndnPacketSize)
{
    throw std::runtime_error("unimplemented method");
}

void ArcModule::dataReceivedX(const std::string &interestName,
                              const std::string &dataName,
                              unsigned int threadId,
			      unsigned int payloadSize,
                              unsigned int ndnPacketSize,
                              int32_t dataNonce,
                              int32_t dGen)
{
    double drd = 0;
    ArcTval tv;
#ifdef ARC_DEBUG_RCVDATA
    std::cout << "ArcModule dataReceived thread_id[" << threadId << "] name[" << interestName << "] size[" << ndnPacketSize << "]" << std::endl;
#endif //ARC_DEBUG_RCVDATA
    getNowTval(&tv);
    sumDataSize_ += ndnPacketSize;

    if (!(consumerPhase_ == ConsumerPhaseFetch || consumerPhase_ == ConsumerPhaseChallenge)) return;
    
    if (threadId == currThreadId_ && currThHist_ != NULL) {
        drd = currThHist_->dataReceivedX(interestName, threadId, payloadSize, ndnPacketSize, dataNonce, dGen);
    } else if (threadId == nextThreadId_ && nextThHist_ != NULL) {
        drd = nextThHist_->dataReceivedX(interestName, threadId, payloadSize, ndnPacketSize, dataNonce, dGen);
    }

#ifdef ARC_EVALUATION
    if (drd != 0)
      //std::cout << "ArcModule " << diffArcTval(&tv, &arcStartTval_) << " Phase " << consumerPhase_ << " DRD " << drd << std::endl;
#endif //ARC_EVALUATION

    return;
}

void ArcModule::updateIndicators(const ArcModule::ArcIndicators& indicators)
{
#ifdef ARC_DEBUG_STATE
    if( (indicators.updateMask_ & BIT_FLAG_7) != 0)
        std::cout << "ArcModule updateIndicators NDNRTC_state:" << consumerPhase_ << "->" << indicators.consumerPhase_ << " state_check_:" << arcStateCheck_ << std::endl;
#endif //ARC_DEBUG_STATE

    if (consumerPhase_ == ConsumerPhaseAdjust && (indicators.updateMask_& BIT_FLAG_7) != 0) {
        if (indicators.consumerPhase_ == ConsumerPhaseFetch) {
	    if(currThHist_ != NULL) {
	        delete currThHist_;
		currThHist_ = NULL;
	    }
            currThHist_ = new ArcHistry;
            arcState_ = arcStateNormal;
	    countChallengePhase_ = 0;
        }
        
    } else if (consumerPhase_ == ConsumerPhaseFetch && (indicators.updateMask_& BIT_FLAG_7) != 0) {
        if (indicators.consumerPhase_ == ConsumerPhaseChallenge) {
            if (arcState_ != onChallengeStarted)
                arcStateCheck_ = false;
	    if(nextThHist_ != NULL) {
	        delete nextThHist_;
		nextThHist_ = NULL;
	    }
            nextThHist_ = new ArcHistry;
            arcState_ = arcStateNormal;

        } else {
	    delete currThHist_;
	    currThHist_ = NULL;
	    arcState_ = arcStateNormal;
	}

    } else if (consumerPhase_ == ConsumerPhaseChallenge && (indicators.updateMask_& BIT_FLAG_7) != 0) {
        if (indicators.consumerPhase_ == ConsumerPhaseFetch) {
            if (arcState_ != onChallengeStopped)
                arcStateCheck_ = false;
	    if(currThHist_ != NULL) {
	        delete currThHist_;
		currThHist_ = NULL;
	    }
	    if(nextThHist_ != NULL) {
	        delete nextThHist_;
		nextThHist_ = NULL;
	    }
            currThHist_ = new ArcHistry;
            arcState_ = arcStateNormal;
	    countChallengePhase_ = 0;

        } else {
	    if(currThHist_ != NULL) {
	        delete currThHist_;
		currThHist_ = NULL;
	    }
	    if(nextThHist_ != NULL) {
	        delete nextThHist_;
		nextThHist_ = NULL;
	    }
	    arcState_ = arcStateNormal;
	    countChallengePhase_ = 0;
	}
    }

    consumerPhase_ = indicators.consumerPhase_;
    return;
}

void ArcModule::reportThreadEvent(const ArcModule::ThreadEvent& event)
{
    if (event ==  ArcModule::ThreadEvent::OldThreadComplete) {
        if (arcState_ != onThreadSwitch)
            arcStateCheck_ = false;
#ifdef ARC_DEBUG_STATE
	std::cout << "ArcModule catches reportThreadEvent:OldThreadComplete state_check:" << arcStateCheck_ << std::endl;
#endif //ARC_DEBUG_STATE
        if (consumerPhase_ == ConsumerPhaseChallenge) {
            arcState_ = onChallengeStopped;
            callback_->onChallengePhaseStopped();
        }
    }
    return;
}


void ArcModule::autoRateControl()
{
    EstResult result_curr, result_next;
    double delta_rate;
    ArcTval tv;
    long interval;
    bool keep_idel_rate = false;
    unsigned int throughput = 0;
    
    getNowTval(&tv);
    interval = diffArcTval(&tv, &arcCallTval_);
#ifdef ARC_DEBUG_TIMER
    std::cout << "ArcModule autoRateControl call interval:" << interval << "[ms]" << std::endl;
#endif //ARC_DEBUG_TIMER
    arcCallTval_ = tv;
    
    if (diffArcTval(&tv, &arcCallThroughputTval_) >= 1000) {
        throughput = (unsigned int)(((double)sumDataSize_ * 8 * 1000) / (diffArcTval(&tv, &arcCallThroughputTval_) * 1024));
	sumDataSize_ = 0;
	arcCallThroughputTval_ = tv;
#ifdef ARC_EVALUATION
        std::cout << "ArcModule " << diffArcTval(&tv, &arcStartTval_) << " Phase " << consumerPhase_ << " Throughput " << throughput << std::endl;
#endif //ARC_EVALUATION
    }

    if (arcState_ == arcStateNormal &&
        (consumerPhase_ == ConsumerPhaseFetch || consumerPhase_ == ConsumerPhaseChallenge)) {

        /* Adaptive Rate Control for Challenge Phase */
        if (consumerPhase_ == ConsumerPhaseChallenge
            && currThHist_ != NULL && nextThHist_ != NULL) {
            result_curr = currThHist_->nwEstimate(interval);
            result_next = nextThHist_->nwEstimate(interval);

	    /* if actual throughtput is not enough than estimate rate, arc module avoid to increase estimate rate */
	    if ( convertPpsToKBps(nextInterestPps_) > (0.8 * nextThHist_->getAvgThroughput ()) ||
		 (convertPpsToKBps(nextInterestPps_) + getBitRateThread(currThreadId_)) > 
		 (0.8 * (currThHist_->getAvgThroughput () + nextThHist_->getAvgThroughput ())))
	        keep_idel_rate = true;

#ifdef ARC_DEBUG_ESTIMATE
	    std::cout << "ArcModule autoRateControl on Challenge Phase thread_id[curr:" << currThreadId_ << " next:" << nextThreadId_ << "] curr_result[" << result_curr << "] next_result [" << result_next << "]" << std::endl;
#endif //ARC_DEBUG_ESTIMATE

            if (result_curr == EstNormal && result_next == EstNormal && keep_idel_rate) {
                nextInterestPps_ += (20 / sqrt (nextInterestPps_));
            } else if (result_curr == EstCongested && result_next == EstCongested) {
                nextInterestPps_ -= (0.5 * sqrt (nextInterestPps_));
            } else if ((result_curr == EstNormal || result_curr == EstUnclear) && result_next == EstCongested) {
                nextInterestPps_ -= (0.5 * sqrt (nextInterestPps_));
            } else if (result_curr == EstCongested && (result_next == EstNormal || result_next == EstUnclear)) {
	        /* there is another option that call onThreadShouldSwitch([lower thread id]) immediately */
                nextInterestPps_ -= (0.5 * sqrt (nextInterestPps_));
            } else if (result_curr == EstCollapse || result_next == EstCollapse) {
	        if (isLowerThread(currThreadId_)) {
		    /* switch thread of lowest bitrate and waiting reportThreadEvent() */
		    callback_->onThreadChallenge(0);
		    currThreadId_ = getLowestThread();
		    arcState_ = onThreadSwitch;
#ifdef ARC_DEBUG_THREAD
		    std::cout << "ArcModule calls onThreadShouldSwitch thread_id[" << currThreadId_ << "] for detecting congestion collapse" << std::endl;
#endif //ARC_DEBUG_THREAD
#ifdef ARC_EVALUATION
		    std::cout << "ArcModule ThreadSwitch thread_id[" << currThreadId_ << "] for detecting congestion collapse" << std::endl;
#endif //ARC_EVALUATION

		    callback_->onThreadShouldSwitch(currThreadId_);
		}
            } else {
                // no process
            }

#ifdef ARC_DEBUG_ESTIMATE
	    std::cout << "ArcModule autoRateControl on Challenge Phase extra_pps[" << nextInterestPps_ << "] extra_bitrate[" << convertPpsToKBps(nextInterestPps_) << "kbp] thread_bitrate[curr:" << getBitRateThread(currThreadId_) << " next:" << getBitRateThread(nextThreadId_) << "]" << std::endl;
#endif //ARC_DEBUG_ESTIMATE
            if ( convertPpsToKBps(nextInterestPps_) >= (STOP_CHALLENGE_RATIO * getBitRateThread(nextThreadId_))) {
	        if (convertPpsToKBps(nextInterestPps_)
                    >= ((1 + SWITCH_HIGHER_SAFTIY_MARGIN) * (getBitRateThread(nextThreadId_) - getBitRateThread(currThreadId_)))) {
                    /* switch thread of higher bitrate and waiting reportThreadEvent() */
                    arcState_ = onThreadSwitch;
                    currThreadId_ = nextThreadId_;
#ifdef ARC_DEBUG_THREAD
                    std::cout << "ArcModule calls onThreadShouldSwitch() thread_id[" << currThreadId_ << "] for switching higher thread" << std::endl;
#endif //ARC_DEBUG_THREAD
#ifdef ARC_EVALUATION
                    std::cout << "ArcModule ThreadSwitch thread_id[" << currThreadId_ << "] for switching higher thread" << std::endl;
#endif //ARC_EVALUATION
                    callback_->onThreadShouldSwitch(currThreadId_);
                } else {
                    delta_rate = convertPpsToKBps(nextInterestPps_) / getBitRateThread(nextThreadId_);
#ifdef ARC_DEBUG_THREAD
                    std::cout << "ArcModule calls onThreadChallenge delta[" << delta_rate << "]" << std::endl;
#endif //ARC_DEBUG_THREAD
                    callback_->onThreadChallenge(delta_rate);
                }
            } else {
	        /* switch thread of lower bitrate and waiting reportThreadEvent() */
	        callback_->onThreadChallenge(0);
                currThreadId_ = getLowerThread(currThreadId_);
                arcState_ = onThreadSwitch;
#ifdef ARC_DEBUG_THREAD
                std::cout << "ArcModule calls onThreadShouldSwitch thread_id[" << currThreadId_ << "] for switching lower thread" << std::endl;
#endif //ARC_DEBUG_THREAD
#ifdef ARC_EVALUATION
                std::cout << "ArcModule ThreadSwitch thread_id[" << currThreadId_ << "] for switching lower thread" << std::endl;
#endif //ARC_EVALUATION
                callback_->onThreadShouldSwitch(currThreadId_);
            }


	/* Adaptive Rate Control for Fetching Phase */
        } else if (consumerPhase_ == ConsumerPhaseFetch && currThHist_ != NULL) {
            result_curr = currThHist_->nwEstimate(interval);
#ifdef ARC_DEBUG_ESTIMATE
	    std::cout << "ArcModule autoRateControl on Fetching Phase thread_id[" << currThreadId_ << "] result[" << result_curr << "]" << std::endl;
#endif //ARC_DEBUG_ESTIMATE
            if (result_curr == EstNormal) {
                ++countChallengePhase_;
                if (countChallengePhase_ >= COUNT_SW_HIGH) {
                    if (isHigherThread(currThreadId_)) {
                        /* start challenge phase and waiting state change notify via updateIndicators() */
                        nextThreadId_ = getHigherThread(currThreadId_);
                        delta_rate = FIRST_CHALLENGE_RATIO * getBitRateThread(nextThreadId_);
                        nextInterestPps_ = convertKBpsToPps(delta_rate);
                        arcState_ = onChallengeStarted;
#ifdef ARC_DEBUG_THREAD
                        std::cout << "ArcModule calls onChallengePhaseStarted thread_id[" << nextThreadId_ << "]" << std::endl;
#endif //ARC_DEBUG_THREAD
                        callback_->onChallengePhaseStarted(nextThreadId_, FIRST_CHALLENGE_RATIO);
                    }
                }
            } else if (result_curr == EstCongested) {
                --countChallengePhase_;
                if (countChallengePhase_ <= COUNT_SW_LOW) {
                    if (isLowerThread(currThreadId_)) {
                        /* switch thread of lower bitrate and waiting reportThreadEvent() */
                        currThreadId_ = getLowerThread(currThreadId_);
                        arcState_ = onThreadSwitch;
#ifdef ARC_DEBUG_THREAD
                        std::cout << "ArcModule calls onThreadShouldSwitch thread_id[" << currThreadId_ << "] for switching lower thread" << std::endl;
#endif //ARC_DEBUG_THREAD
#ifdef ARC_EVALUATION
                        std::cout << "ArcModule ThreadSwitch thread_id[" << currThreadId_ << "] for switching lower thread" << std::endl;
#endif //ARC_EVALUATION
                        callback_->onThreadShouldSwitch(currThreadId_);
                    }
                }
            } else if (result_curr == EstCollapse) {
	        if (isLowerThread(currThreadId_)) {
		  /* switch thread of lowest bitrate and waiting reportThreadEvent() */
		  currThreadId_ = getLowestThread();
		  arcState_ = onThreadSwitch;
#ifdef ARC_DEBUG_THREAD
		  std::cout << "ArcModule calls onThreadShouldSwitch thread_id[" << currThreadId_ << "] for detecting congestion cllapse" << std::endl;
#endif //ARC_DEBUG_THREAD
#ifdef ARC_EVALUATION
		  std::cout << "ArcModule calls ThreadSwitch thread_id[" << currThreadId_ << "] for detecting congestion cllapse" << std::endl;
#endif //ARC_EVALUATION
		  callback_->onThreadShouldSwitch(currThreadId_);
		}
	    }
        }
    }
    
    arcTimer_.expires_from_now(boost::chrono::milliseconds(ARC_INTERVAL));
    arcTimer_.async_wait(boost::bind(&ArcModule::autoRateControl, this));
    return;
}

void ArcModule::getNowTval(ArcTval *qt)
{
    struct tm *tmp;
    struct timeval tv;
    
    if(qt == NULL)return;
    gettimeofday(&tv,NULL);
    tmp=localtime(&tv.tv_sec);
    qt->tim=mktime(tmp);
    qt->msec=tv.tv_usec/1000;
    return;
}


long ArcModule::diffArcTval(const ArcTval* now_t, const ArcTval* prev_t)
{
    long diff_sec;
    long diff_msec;
    diff_sec=difftime(now_t->tim,prev_t->tim);
    diff_msec=now_t->msec-prev_t->msec;
    return diff_sec*1000+diff_msec;
}

double ArcModule::getBitRateThread(const unsigned int threadId)
{
    double oh;
    if (threadId == currThreadId_ && currThHist_ != NULL)
        oh = currThHist_->getHeaderOH ();
    else if (threadId == nextThreadId_ && nextThHist_ != NULL)
        oh = nextThHist_->getHeaderOH ();
    else
        oh = INITIAL_HEADER_OH;

    return (ThreadTable[threadId].bitrate_ * (1.0 + ThreadTable[threadId].parityRatio_) * (1.0 + oh));
}

double ArcModule::convertKBpsToPps(double kbps)
{
    return ((kbps * 1024) / (8 * X_BYTE));
}


double ArcModule::convertPpsToKBps(double pps)
{
    return ((8 * X_BYTE * pps) / 1024);
}


bool ArcModule::isHigherThread(const unsigned int threadId)
{
    ThreadEntry *entry = ThreadTable;
    double curr_rate = ThreadTable[threadId].bitrate_;
    
    for (unsigned int i = 0; i < numThread_; ++i) {
        if (entry->id_ != threadId && entry->bitrate_ > curr_rate)
            return true;
        ++entry;
    }
    return false;
}


bool ArcModule::isLowerThread(const unsigned int threadId)
{
    ThreadEntry *entry = ThreadTable;
    double curr_rate = ThreadTable[threadId].bitrate_;
    
    for (unsigned int i = 0; i < numThread_; ++i) {
        if (entry->id_ != threadId && entry->bitrate_ < curr_rate)
            return true;
        ++entry;
    }
    return false;
}


unsigned int ArcModule::getHigherThread(const unsigned int threadId)
{
    ThreadEntry *entry = ThreadTable;
    double curr_rate = ThreadTable[threadId].bitrate_;
    double next_rate = 0;
    unsigned int target_id = threadId;
    
    for (unsigned int i = 0; i < numThread_; ++i) {
        if (entry->id_ != threadId && entry->bitrate_ > curr_rate) {
            if (next_rate == 0) {
                next_rate = entry->bitrate_;
                target_id = entry->id_;
            } else if (next_rate > 0 && next_rate > entry->bitrate_) {
                next_rate = entry->bitrate_;
                target_id = entry->id_;
            }
        }
        ++entry;
    }
    return target_id;
}


unsigned int ArcModule::getLowerThread(const unsigned int threadId)
{
    ThreadEntry *entry = ThreadTable;
    double curr_rate = ThreadTable[threadId].bitrate_;
    double next_rate = 0;
    unsigned int target_id = threadId;
    
    for (unsigned int i = 0; i < numThread_; ++i) {
        if (entry->id_ != threadId && entry->bitrate_ < curr_rate) {
            if (next_rate == 0) {
                next_rate = entry->bitrate_;
                target_id = entry->id_;
            } else if (next_rate > 0 && next_rate < entry->bitrate_) {
                next_rate = entry->bitrate_;
                target_id = entry->id_;
            }
        }
        ++entry;
    }
    return target_id;
}


unsigned int ArcModule::getLowestThread()
{
    ThreadEntry *entry = ThreadTable;
    double curr_rate = entry->bitrate_;
    unsigned int target_id;
    
    for (unsigned int i = 0; i < numThread_; ++i) {
        if (entry->bitrate_ < curr_rate) {
	    curr_rate = entry->bitrate_;
	    target_id = entry->id_;
        }
        ++entry;
    }
    return target_id;
}


ArcHistry::ArcHistry()
{
    indexSeq_ = lastRcvSeq_ = lastEstSeq_ = 0;
    prevAvgRtt_ = minRtt_ = minRttCandidate_ = 0;
    sumDataSize_ = 0;
    throughtputST_ = throughtputLT_ = 0;
    noRcvCount_ = 0;
    offsetJitter_ = JITTER_OFFSET;
    offsetCollapse_ = COLLAPSE_OFFSET;
    congestionSign_ = false;
}

ArcHistry::~ArcHistry()
{
    
}


void ArcHistry::interestExpressed(const std::string &name,
                                  unsigned int threadId,
				  uint32_t interestNonce)
{
    ArcTval tv;
    getNowTval(&tv);
    ++indexSeq_;
    InterestHistries_.insert (InterestHistry (indexSeq_, name, threadId, tv, interestNonce));
    return;
}


void ArcHistry::interestRetransmit(const std::string &name,
                                   unsigned int threadId)
{
    name_map& nmap = InterestHistries_.get<i_name> ();
    name_map::iterator entry = nmap.find(name);
    // update entry of Interest history
    InterestHistry ih = *entry;
    if (entry == nmap.end ()) return;
    ih.is_retx = true;
    nmap.replace(entry, ih);
    return;
}


double ArcHistry::dataReceivedX(const std::string &name,
                             unsigned int threadId,
			     unsigned int payloadSize,
                             unsigned int ndnPacketSize,
                             uint32_t dataNonce,
                             uint32_t dGen)
{
    ArcTval tv, tv2;
    uint32_t seq, diff_seq;
    double cur_rtt = 0;

    getNowTval(&tv);

    /* counting parameters for caluculate some statstics */
    sumDataSize_ += ndnPacketSize;
    if (listGenDelay_.size () >= LIST_SIZE_GENDLAY)
        listGenDelay_.pop_front ();
    listGenDelay_.push_back (dGen);
    if (listHeaderOH_.size () >= LIST_SIZE_HEADEROH)
        listHeaderOH_.pop_front ();
    listHeaderOH_.push_back ((double)(ndnPacketSize - payloadSize) / ndnPacketSize);

    /* calucurate average of generation delay */
    avgGenDelay_ = calUintListAvg(&listGenDelay_);

    /* finding entry of Interest Histry */
    name_map& nmap = InterestHistries_.get<i_name> ();
    name_map::iterator entry = nmap.find(name);
    if (entry == nmap.end ()) return 0;
    seq = entry->GetSeq ();
    tv2 = entry->GetTxTime ();

    diff_seq = diffSeq (seq, lastRcvSeq_);
    if (diff_seq > 0)
        lastRcvSeq_ = seq;
    
    /* update entry of Interest History */
    InterestHistry ih = *entry;
    ++ih.rx_count;
    if (ih.rx_count == 1) {
        ih.rx_time = tv;
        ih.rtt_prime = diffArcTval(&tv, &tv2);
	//ih.rtt_estimate = diffArcTval(&tv, &tv2) - avgGenDelay_;
	/*
	if (ih.rtt_prime > avgGenDelay_)
	    ih.rtt_estimate = ih.rtt_prime - avgGenDelay_;
	else if (ih.rtt_prime > dGen)
	    ih.rtt_estimate = ih.rtt_prime - dGen;
	*/
	if (ih.rtt_prime > dGen)
	    ih.rtt_estimate = ih.rtt_prime - dGen;
	else
	    ih.rtt_estimate = ih.rtt_prime;
	if (!ih.is_retx) {
	    if (ih.nonce == dataNonce) {
	        ih.is_original = true;
		cur_rtt = ih.rtt_estimate;
	    } else {
	        ih.is_original = false;
		cur_rtt = ih.rtt_prime;
	    }
	}
    }
#ifdef ARC_DEBUG_RCVDATA_DETAIL
    std::cout << "ArcModule dataRecived thread_id[" << ih.GetTid () << "] seq[" << ih.GetSeq () << "] rtt_est[" << ih.GetRttEstimate () << "] rtt_prime[" << ih.GetRttPrime () << "] isOriginal[" << ih.IsOriginal () << "] count[" << ih.GetRxCount () << "] isRetx[" << ih.IsRetx () << "] avgGenDelay[" << avgGenDelay_ << "]" << std::endl;
#endif //ARC_DEBUG_RCVDATA_DETAIL
    nmap.replace(entry, ih);

    /* update minimum rtt */
    if (cur_rtt > 0) {
        long diff_rtt;
        // init param for first data receive
        if (minRtt_ == 0
            && minRttCandidate_ == 0) {
            minRtt_ = minRttCandidate_ = cur_rtt;
            updateMinRttTval_ = tv;
        }
        // update minimum RTT
        if (minRtt_ > cur_rtt)
            minRtt_ = cur_rtt;
        if (minRttCandidate_ > cur_rtt)
            minRttCandidate_ = cur_rtt;
        if (diffArcTval(&tv, &updateMinRttTval_) >  MIN_RTT_EXPIRE) {
            minRtt_ = minRttCandidate_;
            minRttCandidate_ = cur_rtt;
            updateMinRttTval_ = tv;
        }
    }
    return cur_rtt;
}


enum EstResult ArcHistry::nwEstimate(long interval)
{
    uint32_t start_seq = lastEstSeq_ + 1;
    unsigned int rx_count = 0;
    unsigned int loss_count = 0;
    double sum_rtt = 0;
    double avg_rtt = 0;
    double prev_avg_rtt;
    unsigned int tid;
    ArcTval tv, tv2;
    
    getNowTval(&tv);

    /*
    if (diffSeq(lastRcvSeq_, lastEstSeq_) <= 0) {
        ++noRcvCount_;
	if (noRcvCount_ >= ARC_TIMEOUT_COUNT)
	    return EstCollapse;
	else
	    return EstUnclear;
    } else {
        noRcvCount_ = 0;
    }
    */

    seq_map& smap = InterestHistries_.get<i_seq> ();
    
    /* serch entry of Interest Histry for a ARC estimation priod */
    for(uint32_t i = start_seq; i <= lastRcvSeq_; ++i) {
        seq_map::iterator tmp_entry = smap.find(i);
        if (tmp_entry != smap.end ()) {
	    tid = tmp_entry->GetTid ();
            if (!tmp_entry->IsRetx () && tmp_entry->GetRxCount () == 1) {
	        if (tmp_entry->IsOriginal ()) {
		  sum_rtt += tmp_entry->GetRttEstimate ();
		} else {
		  sum_rtt += tmp_entry->GetRttPrime ();
		}
		++rx_count;
            } else {
                ++loss_count;
            }
            lastEstSeq_ = i;
            //delete entry of Interest history
            smap.erase (i);
        }
    }

    /* calcurate average throughput and rtt, and then judge network condition */
    /* if NDN had received data for a ARC estimation priod */
    throughtputST_ = (sumDataSize_ * 8 * 1000) / (interval * 1024);
    listSumDataSize_.push_back (sumDataSize_);
    listArcTval_.push_back (tv);
    tv2 = listArcTval_.front ();
    if (diffArcTval(&tv, &tv2) >= 1000) {
        throughtputLT_ = (unsigned int)((double)calUintListSum(&listSumDataSize_) * 8 * 1000) / (diffArcTval(&tv, &tv2) * 1024);
	listSumDataSize_.pop_front ();
	listArcTval_.pop_front ();
    }
    sumDataSize_ = 0;

    if (rx_count > 0) {
        avg_rtt = sum_rtt / rx_count;

        if (prevAvgRtt_ == 0)
            prevAvgRtt_ = avg_rtt;
	prev_avg_rtt = prevAvgRtt_;
	prevAvgRtt_ = avg_rtt;

#ifdef ARC_DEBUG_ESTIMATE
	std::cout << "ArcModule nwEstimate thread_id[" << tid << "] avg_rtt[" << avg_rtt << "] prev_rtt [" << prev_avg_rtt << "] min_rtt[" << minRtt_ << "] num_of_data[" << rx_count << "] throughput[" << throughtputST_ << "kbps] avg_throughput[" << throughtputLT_ << "kbps]" << std::endl;
#endif //ARC_DEBUG_ESTIMATE

	if (avg_rtt > (minRtt_ + offsetCollapse_))
	    return EstCollapse;

        if (avg_rtt <= (minRtt_ + offsetJitter_) && !congestionSign_) {
            return EstNormal;
        } else if (avg_rtt <= prev_avg_rtt ) {
            if (!congestionSign_) {
                return EstNormal;
            } else {
                congestionSign_ = false;
                return EstUnclear;
            }
        } else {
            if (congestionSign_) {
                return EstCongested;
            } else {
                congestionSign_ = true;
                return EstUnclear;
            }
	}
    }
    return EstUnclear;
}


double ArcHistry::getHeaderOH()
{
    double oh;
    if (listHeaderOH_.size () < (LIST_SIZE_HEADEROH / 4))
        return INITIAL_HEADER_OH;
    oh = calDoubleListAvg(&listHeaderOH_);
    return oh;
}


unsigned int ArcHistry::getAvgThroughput()
{
    return throughtputLT_;
}


unsigned int ArcHistry::getThroughput()
{
    return throughtputST_;
}


void ArcHistry::getNowTval(ArcTval *qt)
{
    struct tm *tmp;
    struct timeval tv;
    
    if(qt == NULL)return;
    gettimeofday(&tv,NULL);
    tmp=localtime(&tv.tv_sec);
    qt->tim=mktime(tmp);
    qt->msec=tv.tv_usec/1000;
    return;
}


long ArcHistry::diffArcTval(const ArcTval* now_t, const ArcTval* prev_t){
    long diff_sec;
    long diff_msec;
    diff_sec=difftime(now_t->tim,prev_t->tim);
    diff_msec=now_t->msec-prev_t->msec;
    return diff_sec*1000+diff_msec;
}


uint32_t ArcHistry::diffSeq (uint32_t a, uint32_t b)
{
    uint32_t diff;
    if (a >= b) {
        diff = a - b;
    }
    else {
        if ((b - a) > (std::numeric_limits<uint32_t>::max () / 2))
            diff = a + (std::numeric_limits<uint32_t>::max () - b);
        else
            diff = a - b;
    }
    return diff;
}



double ArcHistry::calDoubleListAvg (DoubleList *v)
{
    double sum = 0;
    int sample = v->size ();

    if (sample == 0) {
        return 0;
    }
    std::list<double>::iterator TLIterator = v->begin ();
    while (TLIterator != v->end ()) {
        sum += *TLIterator;
	TLIterator++;
    }
    return (sum / sample);
}


unsigned int  ArcHistry::calUintListAvg (UintList *v)
{
    unsigned sum = 0;
    int sample = v->size ();

    if (sample == 0) {
        return 0;
    }
    std::list<unsigned int>::iterator TLIterator = v->begin ();
    while (TLIterator != v->end ()) {
        sum += *TLIterator;
	TLIterator++;
    }
    return (sum / sample);
}


unsigned int  ArcHistry::calUintListSum (UintList *v)
{
    unsigned sum = 0;
    int sample = v->size ();

    if (sample == 0) {
        return 0;
    }
    std::list<unsigned int>::iterator TLIterator = v->begin ();
    while (TLIterator != v->end ()) {
        sum += *TLIterator;
	TLIterator++;
    }
    return sum;
}
