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
    getNowTval(&arcCallTval_);
    
    for (auto itr = mediaThreads.begin(); itr != mediaThreads.end(); ++itr) {
        ThreadTable[numThread_] = *itr;
#ifdef ARC_DEBUG
        std::cout << "ArcModule initialize() thread_id[" << ThreadTable[numThread_].id_ << "] bitrate[" <<  ThreadTable[numThread_].bitrate_ << "] parity[" << ThreadTable[numThread_].parityRatio_ << "]" << std::endl;
#endif //ARC_DEBUG
        ++numThread_;
    }
    return 0;
}

void ArcModule::interestExpressed(const std::string &name,
                                  unsigned int threadId,
				  uint32_t interestNonce)
{
#ifdef ARC_DEBUG
    std::cout << "ArcModule interestExpressed() thread_id[" << threadId << "] name[" << name << "]" << std::endl;
#endif //ARC_DEBUG
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
#ifdef ARC_DEBUG
    std::cout << "ArcModule interestRetransmit() thread_id[" << threadId << "] name[" << name << "]" << std::endl;
#endif //ARC_DEBUG
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
                              unsigned int ndnPacketSize,
                              int32_t dataNonce,
                              int32_t dGen)
{
#ifdef ARC_DEBUG
    std::cout << "ArcModule dataReceived() thread_id[" << threadId << "] name[" << interestName << "] size[" << ndnPacketSize << "]" << std::endl;
#endif //ARC_DEBUG
    if (!(consumerPhase_ == ConsumerPhaseFetch || consumerPhase_ == ConsumerPhaseChallenge)) return;
    
    if (threadId == currThreadId_ && currThHist_ != NULL) {
        currThHist_->dataReceived(interestName, threadId, ndnPacketSize, dataNonce, dGen);
    } else if (threadId == nextThreadId_ && nextThHist_ != NULL) {
        nextThHist_->dataReceived(interestName, threadId, ndnPacketSize, dataNonce, dGen);
    }
    return;
}

void ArcModule::updateIndicators(const ArcModule::ArcIndicators& indicators)
{
    if (consumerPhase_ == ConsumerPhaseAdjust) {
        if (indicators.consumerPhase_ == ConsumerPhaseFetch) {
            std::cout << "ArcModule updateIndicators() NDNRTC state:Adjust->Fetch state_check:" << arcStateCheck_ << std::endl;
#ifdef ARC_DEBUG
            std::cout << "ArcModule updateIndicators() NDNRTC state:Adjust->Fetch state_check:" << arcStateCheck_ << std::endl;
#endif //ARC_DEBUG
            currThHist_ = new ArcHistry;
            arcState_ = arcStateNormal;
	    countChallengePhase_ = 0;
        }
        
    } else if (consumerPhase_ == ConsumerPhaseFetch) {
        if (indicators.consumerPhase_ == ConsumerPhaseAdjust) {
            std::cout << "ArcModule updateIndicators() NDNRTC state:Fetch->Adjust state_check:" << arcStateCheck_ << std::endl;
#ifdef ARC_DEBUG
            std::cout << "ArcModule updateIndicators() NDNRTC state:Fetch->Adjust state_check:" << arcStateCheck_ << std::endl;
#endif //ARC_DEBUG
            delete currThHist_;
            currThHist_ = NULL;
            arcState_ = arcStateNormal;
        }
        if (indicators.consumerPhase_ == ConsumerPhaseChallenge) {
            if (arcState_ != onChallengeStarted)
                arcStateCheck_ = false;
            std::cout << "ArcModule updateIndicators() NDNRTC state:Fetch->Challenge state_check:" << arcStateCheck_ << std::endl;
#ifdef ARC_DEBUG
            std::cout << "ArcModule updateIndicators() NDNRTC state:Fetch->Challenge state_check:" << arcStateCheck_ << std::endl;
#endif //ARC_DEBUG
            nextThHist_ = new ArcHistry;
            arcState_ = arcStateNormal;
        }
        
    } else if (consumerPhase_ == ConsumerPhaseChallenge) {
        if (indicators.consumerPhase_ == ConsumerPhaseFetch) {
            if (arcState_ != onChallengeStopped)
                arcStateCheck_ = false;
            std::cout << "ArcModule updateIndicators() NDNRTC state:Challenge->Fetch state_check:" << arcStateCheck_ << std::endl;
#ifdef ARC_DEBUG
            std::cout << "ArcModule updateIndicators() NDNRTC state:Challenge->Fetch state_check:" << arcStateCheck_ << std::endl;
#endif //ARC_DEBUG
            delete currThHist_;
            currThHist_ = new ArcHistry;
            nextThHist_ = new ArcHistry;
            arcState_ = arcStateNormal;
	    countChallengePhase_ = 0;
        }
    }
    
    consumerPhase_ = indicators.consumerPhase_;
    return;
}

void ArcModule::reportThreadEvent(const ArcModule::ThreadEvent& event)
{
#ifdef ARC_DEBUG
    std::cout << "ArcModule reportThreadEvent() state_check:" << arcStateCheck_ << std::endl;
#endif //ARC_DEBUG
    if (event ==  ArcModule::ThreadEvent::OldThreadComplete) {
        if (arcState_ != onThreadSwitch)
            arcStateCheck_ = false;
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
    
    getNowTval(&tv);
#ifdef ARC_DEBUG
    std::cout << "ArcModule autoRateControl() interval:" << diffArcTval(&tv, &arcCallTval_) << "[ms]" << std::endl;
#endif //ARC_DEBUG
    arcCallTval_ = tv;
    
    if (arcState_ == arcStateNormal &&
        (consumerPhase_ == ConsumerPhaseFetch || consumerPhase_ == ConsumerPhaseChallenge)) {
        
        if (consumerPhase_ == ConsumerPhaseChallenge
            && currThHist_ != NULL && nextThHist_ != NULL) {
            result_curr = currThHist_->nwEstimate();
            result_next = nextThHist_->nwEstimate();

            if (result_curr == EstNormal && result_next == EstNormal) {
                nextInterestPps_ += (20 / sqrt (nextInterestPps_));
            } else if (result_curr == EstCongested && result_next == EstCongested) {
                nextInterestPps_ -= (0.5 * sqrt (nextInterestPps_));
            } else if (result_curr == EstNormal && result_next == EstCongested) {
                nextInterestPps_ -= (0.5 * sqrt (nextInterestPps_));
            } else if (result_curr == EstCongested && result_next == EstNormal) {
                nextInterestPps_ -= (0.5 * sqrt (nextInterestPps_));
                // there is another option that call onThreadShouldSwitch([lower thread id]) immediately
            } else if (result_curr == EstCollapse || result_next == EstCollapse) {
	      if (isLowerThread(currThreadId_)) {
		/* switch thread of lower bitrate and waiting reportThreadEvent() */
		callback_->onThreadChallenge(0);
		currThreadId_ = getLowerThread(currThreadId_);
		arcState_ = onThreadSwitch;
#ifdef ARC_DEBUG
		std::cout << "ArcModule autoRateControl() call onThreadShouldSwitch() thread_id[" << currThreadId_ << "]" << std::endl;
#endif //ARC_DEBUG
		callback_->onThreadShouldSwitch(currThreadId_);
	      }
            } else {
                // no process
            }

	    //std::cout << << std::endl;

            if (nextInterestPps_ > 0) {
                if (convertPpsToKBps(nextInterestPps_)
                    >= getBitRateThread(nextThreadId_) - getBitRateThread(currThreadId_)) {
                    /* switch thread of higher bitrate and waiting reportThreadEvent() */
                    arcState_ = onThreadSwitch;
                    currThreadId_ = nextThreadId_;
#ifdef ARC_DEBUG
                    std::cout << "ArcModule autoRateControl() call onThreadShouldSwitch() thread_id[" << currThreadId_ << "]" << std::endl;
#endif //ARC_DEBUG
                    callback_->onThreadShouldSwitch(currThreadId_);
                } else {
                    delta_rate = convertPpsToKBps(nextInterestPps_) / getBitRateThread(nextThreadId_);
#ifdef ARC_DEBUG
                    std::cout << "ArcModule autoRateControl() call onThreadChallenge() delta[" << delta_rate << "]" << std::endl;
#endif //ARC_DEBUG
                    callback_->onThreadChallenge(delta_rate);
                }
            } else if (nextInterestPps_ <= 0) {
                /* switch thread of lower bitrate and waiting reportThreadEvent() */
                //callback_->onThreadChallenge(0);
                currThreadId_ = getLowerThread(currThreadId_);
                arcState_ = onThreadSwitch;
#ifdef ARC_DEBUG
                std::cout << "ArcModule autoRateControl() call onThreadShouldSwitch() thread_id[" << currThreadId_ << "]" << std::endl;
#endif //ARC_DEBUG
                callback_->onThreadShouldSwitch(currThreadId_);
            }
            
        } else if (consumerPhase_ == ConsumerPhaseFetch && currThHist_ != NULL) {
            result_curr = currThHist_->nwEstimate();
            if (result_curr == EstNormal) {
                ++countChallengePhase_;
                if (countChallengePhase_ >= COUNT_SW_HIGH) {
                    if (isHigherThread(currThreadId_)) {
                        /* start challenge phase and waiting state change notify via updateIndicators() */
                        nextThreadId_ = getHigherThread(currThreadId_);
                        delta_rate = FIRST_CHALLENGE_RATIO * getBitRateThread(nextThreadId_);
                        nextInterestPps_ = FIRST_CHALLENGE_RATIO * convertKBpsToPps(delta_rate);
                        arcState_ = onChallengeStarted;
#ifdef ARC_DEBUG
                        std::cout << "ArcModule autoRateControl() call onChallengePhaseStarted() thread_id[" << nextThreadId_ << "]" << std::endl;
#endif //ARC_DEBUG
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
#ifdef ARC_DEBUG
                        std::cout << "ArcModule autoRateControl() call onThreadShouldSwitch() thread_id[" << currThreadId_ << "]" << std::endl;
#endif //ARC_DEBUG
                        callback_->onThreadShouldSwitch(currThreadId_);
                    }
                }
            } else if (result_curr == EstCollapse) {
	        if (isLowerThread(currThreadId_)) {
		  /* switch thread of lower bitrate and waiting reportThreadEvent() */
		  currThreadId_ = getLowerThread(currThreadId_);
		  arcState_ = onThreadSwitch;
#ifdef ARC_DEBUG
		  std::cout << "ArcModule autoRateControl() call onThreadShouldSwitch() thread_id[" << c\
		    urrThreadId_ << "]" << std::endl;
#endif //ARC_DEBUG                                                                                             
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
    //printf("%04d/%02d/%02d %02d:%02d:%02d:%3d\n",
    //      tmp->tm_year + 1900, tmp->tm_mon + 1,
    //      tmp->tm_mday, tmp->tm_hour,
    //      tmp->tm_min, tmp->tm_sec,
    //      tv.tv_usec/1000);
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
    ThreadEntry *entry = ThreadTable;
    for (unsigned int i = 0; i < numThread_; ++i) {
        if (entry->id_ == threadId)
            return (entry->bitrate_ * (1.0 + entry->parityRatio_));
        ++entry;
    }
    return 0;
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
    double curr_rate = getBitRateThread(threadId);
    
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
    double curr_rate = getBitRateThread(threadId);
    
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
    double curr_rate = getBitRateThread(threadId);
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
    double curr_rate = getBitRateThread(threadId);
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


ArcHistry::ArcHistry()
{
    indexSeq_ = lastRcvSeq_ = lastEstSeq_ = 0;
    prevAvgRtt_ = minRtt_ = minRttCandidate_ = 0;
    avgDataSize_ = 0;
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


void ArcHistry::dataReceived(const std::string &name,
                             unsigned int threadId,
                             unsigned int ndnPacketSize,
                             uint32_t dataNonce,
                             uint32_t dGen)
{
    ArcTval tv, tv2;
    uint32_t seq, diff_seq;
    double cur_rtt;

    getNowTval(&tv);
    
    if (avgDataSize_ == 0)
        avgDataSize_ = ndnPacketSize;
    avgDataSize_ = 0.9 * avgDataSize_ + 0.1 * ndnPacketSize;
    
    name_map& nmap = InterestHistries_.get<i_name> ();
    name_map::iterator entry = nmap.find(name);
    if (entry == nmap.end ()) return;
    seq = entry->GetSeq ();
    diff_seq = diffSeq (seq, lastRcvSeq_);
    if (diff_seq > 0)
        lastRcvSeq_ = seq;
    
    tv2 = entry->GetTxTime ();
    if (entry->GetNonce() == dataNonce)
        cur_rtt = diffArcTval(&tv, &tv2) - dGen;
    else
        cur_rtt = diffArcTval(&tv, &tv2);

    // update entry of Interest history
    InterestHistry ih = *entry;
    ++ih.rx_count;
    if (ih.rx_count == 1) {
        ih.rx_time = tv;
        ih.rtt_prime = diffArcTval(&tv, &tv2);
	ih.rtt_estimate = diffArcTval(&tv, &tv2) - dGen;
	if (ih.nonce == dataNonce)
	  ih.is_original = true;
	else
	  ih.is_original = false;
    }
    nmap.replace(entry, ih);

    // update minimum rtt
    if (cur_rtt > 0) {
        long diff_rtt;
        // init param for first data receive
        if (minRtt_ == 0
            && minRttCandidate_ == 0) {
            minRtt_ = minRttCandidate_ = cur_rtt;
            updateMinRttTval_ = tv;
        }
        // update minimum RTT
        if (minRttCandidate_ > cur_rtt)
            minRttCandidate_ = cur_rtt;
        if (diffArcTval(&tv, &updateMinRttTval_) >  MIN_RTT_EXPIRE) {
            minRtt_ = minRttCandidate_;
            minRttCandidate_ = cur_rtt;
            updateMinRttTval_ = tv;
        }
    }
    return;
}


enum EstResult ArcHistry::nwEstimate()
{
    uint32_t start_seq = lastEstSeq_ + 1;
    unsigned int rx_count = 0;
    unsigned int loss_count = 0;
    double sum_rtt = 0;
    double avg_rtt = 0;
    double prev_avg_rtt;
    
    if (diffSeq(lastRcvSeq_, lastEstSeq_) <= 0) return EstUnclear;
    
    seq_map& smap = InterestHistries_.get<i_seq> ();
    
    for(uint32_t i = start_seq; i <= lastRcvSeq_; ++i) {
        seq_map::iterator tmp_entry = smap.find(i);
        if (tmp_entry != smap.end ()) {
            if (!tmp_entry->IsRetx ()) {
	        if (tmp_entry->IsOriginal ()) {
		  sum_rtt += tmp_entry->GetRttEstimate ();
		  ++rx_count;
		} else {
		  sum_rtt += tmp_entry->GetRttPrime ();
		  ++rx_count;
		}
            } else {
                ++loss_count;
            }
            lastEstSeq_ = i;
            //delete entry of Interest history
            smap.erase (i);
        }
    }
    
    if (rx_count > 0) {
        avg_rtt = sum_rtt / rx_count;

	//std::cout << "avg_rtt " << avg_rtt << " min_rtt " << minRtt_ << std::endl;

        if (prevAvgRtt_ == 0)
            prevAvgRtt_ = avg_rtt;
	prev_avg_rtt = prevAvgRtt_;
	prevAvgRtt_ = avg_rtt;

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


void ArcHistry::getNowTval(ArcTval *qt)
{
    struct tm *tmp;
    struct timeval tv;
    
    if(qt == NULL)return;
    gettimeofday(&tv,NULL);
    tmp=localtime(&tv.tv_sec);
    qt->tim=mktime(tmp);
    qt->msec=tv.tv_usec/1000;
    //printf("%04d/%02d/%02d %02d:%02d:%02d:%3d\n",
    //      tmp->tm_year + 1900, tmp->tm_mon + 1,
    //      tmp->tm_mday, tmp->tm_hour,
    //      tmp->tm_min, tmp->tm_sec,
    //      tv.tv_usec/1000);
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
