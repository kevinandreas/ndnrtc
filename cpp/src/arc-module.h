//
//  arc-module.h
//  ndnrtc
//
//  Created by Takahiro Yoneda on 11/13/15.
//  Copyright 2013-2015 Regents of the University of California
//

#define THREAD_MAX 10
#define COUNT_SW_HIGH 10
#define COUNT_SW_LOW 10
#define X_BYTE 1024
#define JITTER_OFFSET 10
//#define ARC_DEBUG
#define ARC_INTERVAL 50
#define MIN_RTT_EXPIRE 30000

#include <time.h>
#include <sys/time.h>

#include <boost/asio/steady_timer.hpp>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/tag.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>


#include "arc-interface.h"

namespace ndnrtc {
    
    class ArcHistry;
    
    enum EstResult {
        EstNormal,
        EstUnclear,
        EstCongested,
    };
    
    struct ArcTval {
        time_t tim;
        long msec;
    };
    
    /**
     * A base class for ARC module class which uses boost::asio::io_service
     * for scheduling callbacks and timers
     */
    class ArcModuleIoBase : public IRateAdaptationModule
    {
    public:
        ArcModuleIoBase(boost::asio::io_service& ioService):io_(ioService){}
        virtual ~ArcModuleIoBase(){}
        
    protected:
        boost::asio::io_service& io_;
    };
    
    /**
     * ARC module
     */
    class ArcModule : public ArcModuleIoBase
    {
    public:
        static const ArcIndicators ZeroIndicators;
        
        ArcModule(boost::asio::io_service& ioService):ArcModuleIoBase(ioService), arcTimer_(io_){}
        ~ArcModule() {}
        
        int initialize(IRateAdaptationModuleCallback* const callback,
                       const CodecMode& codecMode,
                       std::vector<ThreadEntry> mediaThreads);
        void interestExpressed(const std::string &name,
                               unsigned int threadId,
			       uint32_t interestNonce);
        void interestRetransmit(const std::string &name,
                                unsigned int threadId);
        void dataReceivedX(const std::string &interestName,
			   const std::string &dataName,
			   unsigned int threadId,
			   unsigned int ndnPacketSize,
			   uint32_t dataNonce,
			   uint32_t dGen);
        void updateIndicators(const ArcIndicators& indicators);
        void reportThreadEvent(const ThreadEvent& event);
        void autoRateControl();
        
    private:
        IRateAdaptationModuleCallback* callback_;
        enum ArcState {
            arcStateNormal,
            onChallengeStarted,
            onThreadSwitch,
            onChallengeStopped,
        };
        
        boost::asio::steady_timer arcTimer_;
        
        ThreadEntry ThreadTable[THREAD_MAX];
        unsigned int numThread_;
        unsigned int currThreadId_;
        unsigned int nextThreadId_;
        enum CodecMode codecMode_;
        ConsumerPhase consumerPhase_;
        enum ArcState arcState_;
        bool arcStateCheck_; /* true:normal false:state error */
        ArcTval arcCallTval_;
        
        double nextInterestPps_;
        ArcHistry *currThHist_, *nextThHist_;
        int countChallengePhase_;
        
        void getNowTval(ArcTval *qt);
        long diffArcTval(const ArcTval* now_t, const ArcTval* prev_t);
        double getBitRateThread(const unsigned int threadId);
        double convertBpsToPps(double bps);
        double convertPpsToBps(double pps);
        bool isHigherThread(const unsigned int threadId);
        bool isLowerThread(const unsigned int threadId);
        unsigned int getHigherThread(const unsigned int threadId);
        unsigned int getLowerThread(const unsigned int threadId);
    };
    
    class ArcHistry
    {
    public:
        ArcHistry ();
        ~ArcHistry ();
        void interestExpressed(const std::string &name,
                               unsigned int threadId,
			       uint32_t interestNonce);
        void interestRetransmit(const std::string &name,
                                unsigned int threadId);
        void dataReceived(const std::string &name,
                          unsigned int threadId,
                          unsigned int ndnPacketSize,
			  uint32_t dataNonce,
			  uint32_t dGen);
        enum EstResult nwEstimate();
        
    private:
        uint32_t indexSeq_, lastRcvSeq_, lastEstSeq_;
        double prevAvgRtt_, minRtt_, minRttCandidate_;
        double avgDataSize_;
        long offsetJitter_;
        ArcTval updateMinRttTval_;
        bool congestionSign_;
        
        uint32_t diffSeq (uint32_t a, uint32_t b);
        void getNowTval(ArcTval *qt);
        long diffArcTval(const ArcTval* now_t, const ArcTval* prev_t);
        
        struct InterestHistry {
            InterestHistry (const uint32_t _seq, const std::string& _name,
                            const unsigned int _tid, const ArcTval _txtime,
			    const uint32_t _nonce) :
	    seq (_seq), name (_name), tid (_tid), tx_time (_txtime), nonce (_nonce),
	    rx_count (0), is_retx (false), rtt_prime (0), rtt_estimate (0),
	    is_original (false) { };
            uint32_t seq;
            std::string name;
            unsigned int tid;
            ArcTval tx_time;
            ArcTval rx_time;
            bool is_retx;
            long rtt_prime;
	    long rtt_estimate;
	    uint32_t nonce;
	    bool is_original;
            unsigned int rx_count;
            uint32_t GetSeq () const { return seq; }
            std::string GetName () const { return name; }
            unsigned int GetTid () const { return tid; }
            ArcTval GetTxTime () const { return tx_time; }
            ArcTval GetRxTime () const { return rx_time; }
	    bool IsRetx () const { return is_retx; }
            long GetRttPrime () const { return rtt_prime; }
            long GetRttEstimate () const { return rtt_estimate; }
	    uint32_t GetNonce () const { return nonce; }
	    bool IsOriginal () const { return is_original; }
            unsigned int GetRxCount () const { return rx_count; }
        };
        
        class i_seq { };
        class i_name { };
        
        typedef boost::multi_index::multi_index_container<
        InterestHistry,
        boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<
        boost::multi_index::tag<i_seq>,
        boost::multi_index::member<InterestHistry, uint32_t, &InterestHistry::seq>
        >,
        boost::multi_index::ordered_unique<
        boost::multi_index::tag<i_name>,
        boost::multi_index::member<InterestHistry, std::string, &InterestHistry::name>
        >
        >
        > InterestHistryContainer;
        typedef InterestHistryContainer::index<i_name>::type name_map;
        typedef InterestHistryContainer::index<i_seq>::type seq_map;
        InterestHistryContainer InterestHistries_;
    };
};
