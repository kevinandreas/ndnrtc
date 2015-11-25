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
    arcTimer_.expires_from_now(boost::chrono::milliseconds(10));
    arcTimer_.async_wait(boost::bind(&ArcModule::autoRateControl, this));
    //arcTimer_.async_wait(&ArcModule::autoRateControl);
    //arcTimer_.async_wait([this, 10, autoRateControl]);
    
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
                                  unsigned int threadId)
{
#ifdef ARC_DEBUG
    std::cout << "ArcModule interestExpressed() thread_id[" << threadId << "] name[" << name << "]" << std::endl;
#endif //ARC_DEBUG
    if (consumerPhase_ == ConsumerPhaseAdjust)
        currThreadId_ = threadId;
    if (!(consumerPhase_ == ConsumerPhaseFetch || consumerPhase_ == ConsumerPhaseChallenge)) return;
    
    if (threadId == currThreadId_ && currThHist_ != NULL) {
        currThHist_->interestExpressed(name, threadId);
    } else if (threadId == nextThreadId_ && nextThHist_ != NULL) {
        nextThHist_->interestExpressed(name, threadId);
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
#ifdef ARC_DEBUG
    std::cout << "ArcModule dataReceived() thread_id[" << threadId << "] name[" << interestName << "] size[" << ndnPacketSize << "]" << std::endl;
#endif //ARC_DEBUG
    if (!(consumerPhase_ == ConsumerPhaseFetch || consumerPhase_ == ConsumerPhaseChallenge)) return;
    
    if (threadId == currThreadId_ && currThHist_ != NULL) {
        currThHist_->dataReceived(interestName, threadId, ndnPacketSize);
    } else if (threadId == nextThreadId_ && nextThHist_ != NULL) {
        nextThHist_->dataReceived(interestName, threadId, ndnPacketSize);
    }
    return;
}

void ArcModule::updateIndicators(const ArcModule::ArcIndicators& indicators)
{
    if (consumerPhase_ == ConsumerPhaseAdjust) {
        if (indicators.consumerPhase_ == ConsumerPhaseFetch) {
#ifdef ARC_DEBUG
            std::cout << "ArcModule updateIndicators() NDNRTC state:Adjust->Fetch state_check:" << arcStateCheck_ << std::endl;
#endif //ARC_DEBUG
            currThHist_ = new ArcHistry;
            arcState_ = arcStateNormal;
        }
        
    } else if (consumerPhase_ == ConsumerPhaseFetch) {
        countChallengePhase_ = 0;
        if (indicators.consumerPhase_ == ConsumerPhaseAdjust) {
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
#ifdef ARC_DEBUG
            std::cout << "ArcModule updateIndicators() NDNRTC state:Challenge->Fetch state_check:" << arcStateCheck_ << std::endl;
#endif //ARC_DEBUG
            delete currThHist_;
            currThHist_ = new ArcHistry;
            nextThHist_ = new ArcHistry;
            arcState_ = arcStateNormal;
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
//void ArcModule::autoRateControl(const boost::system::error_code &event)
{
    EstResult result_curr, result_next;
    double delta_rate;
    ArcTval tv;
    
    /*
     if (event == boost::asio::error::operation_aborted) {
     std::cout << "ArcModule autoRateControl() timer error" << std::endl;
     return;
     }
     */
    
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
            } else {
                // no process
            }
            
            if (nextInterestPps_ > 0) {
                if (convertPpsToBps(nextInterestPps_)
                    >= getBitRateThread(nextThreadId_) - getBitRateThread(currThreadId_)) {
                    /* switch thread of higher bitrate and waiting reportThreadEvent() */
                    arcState_ = onThreadSwitch;
                    currThreadId_ = nextThreadId_;
#ifdef ARC_DEBUG
                    std::cout << "ArcModule autoRateControl() call onThreadShouldSwitch() thread_id[" << currThreadId_ << "]" << std::endl;
#endif //ARC_DEBUG
                    callback_->onThreadShouldSwitch(currThreadId_);
                } else {
                    delta_rate = convertPpsToBps(nextInterestPps_) / getBitRateThread(nextThreadId_);
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
                        delta_rate = 0.05 * getBitRateThread(nextThreadId_);
                        nextInterestPps_ = 0.05 * convertBpsToPps(delta_rate);
                        arcState_ = onChallengeStarted;
#ifdef ARC_DEBUG
                        std::cout << "ArcModule autoRateControl() call onChallengePhaseStarted() thread_id[" << nextThreadId_ << "]" << std::endl;
#endif //ARC_DEBUG
                        callback_->onChallengePhaseStarted(nextThreadId_, nextInterestPps_);
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
            }
        }
    }
    
    arcTimer_.expires_from_now(boost::chrono::milliseconds(10));
    //arcTimer_.async_wait(&ArcModule::autoRateControl);
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

double ArcModule::convertBpsToPps(double bps)
{
    return (bps / (8 * X_BYTE));
}


double ArcModule::convertPpsToBps(double pps)
{
    return (8 * X_BYTE * pps);
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
    prevAvgRtt_ = minRtt_ = minRttCandidate_ = std::numeric_limits<long>::max ();
    avgDataSize_ = 0;
    offsetJitter_ = JITTER_OFFSET;
    congestionSign_ = false;
}

ArcHistry::~ArcHistry()
{
    
}


void ArcHistry::interestExpressed(const std::string &name,
                                  unsigned int threadId)
{
    ArcTval tv;
    getNowTval(&tv);
    ++indexSeq_;
    InterestHistries_.insert (InterestHistry (indexSeq_, name, threadId, tv));
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
    ih.retx_flag = true;
    nmap.replace(entry, ih);
    return;
}


void ArcHistry::dataReceived(const std::string &name,
                             unsigned int threadId,
                             unsigned int ndnPacketSize)
{
    ArcTval tv, tv2;
    uint32_t seq, diff_seq;
    long cur_rtt = 0;
    
    getNowTval(&tv);
    
    if (avgDataSize_ == 0)
        avgDataSize_ = ndnPacketSize;
    avgDataSize_ = 0.9 * avgDataSize_ + 0.1 * ndnPacketSize;
    
    name_map& nmap = InterestHistries_.get<i_name> ();
    name_map::iterator entry = nmap.find(name);
    if (entry == nmap.end ()) return;
    seq = entry->GetSeq ();
    tv2 = entry->GetTxTime ();
    diff_seq = diffSeq (seq, lastRcvSeq_);
    if (diff_seq > 0)
        lastRcvSeq_ = seq;
    cur_rtt = diffArcTval(&tv, &tv2);
    
    // update entry of Interest history
    InterestHistry ih = *entry;
    ++ih.rx_count;
    if (ih.rx_count == 1) {
        ih.rx_time = tv;
        ih.rtt = cur_rtt;
    }
    nmap.replace(entry, ih);
    
    // update minimum rtt
    if (cur_rtt > 0) {
        long diff_rtt;
        // init param for first data receive
        if (minRtt_ == std::numeric_limits<long>::max ()
            && minRttCandidate_ == std::numeric_limits<long>::max ()) {
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
    long sum_rtt = 0;
    long avg_rtt = 0;
    
    if (diffSeq(lastRcvSeq_, lastEstSeq_) <= 0) return EstUnclear;
    
    seq_map& smap = InterestHistries_.get<i_seq> ();
    
    for(uint32_t i = start_seq; i <= lastRcvSeq_; ++i) {
        seq_map::iterator tmp_entry = smap.find(i);
        if (tmp_entry != smap.end ()) {
            if (tmp_entry->GetRxCount () > 0) {
                sum_rtt += tmp_entry->GetRtt ();
                ++rx_count;
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
        if (prevAvgRtt_ == std::numeric_limits<long>::max ())
            prevAvgRtt_ = avg_rtt;
        if (avg_rtt <= (minRtt_ + offsetJitter_) && loss_count == 0 && !congestionSign_) {
            return EstNormal;
        } else if (avg_rtt <= prevAvgRtt_ && loss_count == 0) {
            if (!congestionSign_) {
                return EstNormal;
            } else {
                congestionSign_ = false;
                return EstUnclear;
            }
        } else if (avg_rtt > prevAvgRtt_ && loss_count == 0) {
            if (congestionSign_) {
                return EstCongested;
            } else {
                congestionSign_ = true;
                return EstUnclear;
            }
        } else {
            // if rx_count > 0 && loss_count > 0
            congestionSign_ = true;
            return EstCongested;
        }
        prevAvgRtt_ = avg_rtt;
    } else if (loss_count > 0) {
        congestionSign_ = true;
        return EstCongested;
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
