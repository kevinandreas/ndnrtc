//
//  arc-interface.h
//  ndnrtc
//
//  Created by Peter Gusev on 11/13/15.
//  Copyright 2013-2015 Regents of the University of California
//
//  Authors:  Peter Gusev, Ryota Ohnishi, Takahiro Yoneda
//

#ifndef arc_interface_h
#define arc_interface_h

#include <string>
#include <stdint.h>
#include <boost/asio.hpp>

namespace ndnrtc {

    class IRateAdaptationModuleCallback;

    /**
     * This abstract class defines an interface for rate adaptation decision
     * module which provides several functions:
     *  - detects network congestions;
     *  - recommends video bitrate for the current network conditions;
     *  - recommends interest issuing rate for the current network conditions.
     */
    class IRateAdaptationModule {
    	public:
     	/**
	      * Codec modes
	      */
    	enum class CodecMode {
        	CodecModeNormal, // normal codec mode
        	CodecModeSVC     // scalable video coding (SVC) - not supproted
	    };
     
     	typedef struct _ThreadEntry {
        	int id_;
         	double bitrate_;
         	double parityRatio_;
     	} ThreadEntry;
     
     	typedef enum _ConsumerPhase {  
     		ConsumerPhaseInactive,
     		ConsumerPhaseWaitInitial,
     		ConsumerPhaseChasing,
     		ConsumerPhaseAdjust,
     		ConsumerPhaseFetch,
     		ConsumerPhaseChallenge
     	} ConsumerPhase;

     	typedef struct _ArcIndicators {
     		unsigned int updateMask_;	// if bit is set, corresponding value has new data
     									// see comments below for which value corresponds
     									// to which bit number
     		double Darr_;				// bit 0 in updateMask
     		double producerRate_;		// bit 1 in updateMask
     		double rttPrime_;			// bit 2 in updateMask
     		double rttRealEstimated_;	// bit 3 in updateMask
     		double bufferTargetSize_;	// bit 4 in updateMask
     		double bufferPlayableSize_;	// bit 5 in updateMask
     		double bufferReservedSize_; // bit 6 in updateMask
     		ConsumerPhase consumerPhase_; // bit 7 in updateMask
     	} ArcIndicators;
     
     	enum class ThreadEvent {
     		InterestsSwitched,	// when consumer started to issue interests towards new thread
     		NewThreadStarted,	// when data for new thread begins to arrive
     		OldThreadComplete	// when last data for old thread has arrived
     	};
     	
     	/**
     	 * ARC module initialization
     	 * @param callback - 	ARC module callback - represents NDN-RTC for ARC module
     	 * @param codecMode - 	Encoder mode (@see CodecMode)
     	 * @param mediaThreads -Media threads available for consumer to fetch
     	 * @return 0 if initialization was successfull, non-zero otherwise
     	 */
        virtual int initialize(IRateAdaptationModuleCallback* const callback,
        						const CodecMode& codecMode,
        						std::vector<ThreadEntry> mediaThreads) = 0;
        /**
         * Called by NDN-RTC every time new Interest is issued
         * @param name -	Interest name
         * @param threadId -ID of the media thread issued Interest belongs to
         * @param interestNonce - Interest nonce value
         */
        virtual void interestExpressed(const std::string &name,
                                       unsigned int threadId,
                                       uint32_t interestNonce) = 0;

        /**
         * Called by NDN-RTC every time Interest retransmission happened
         * @param name -	Interest name
         * @param threadId -ID of the media thread issued Interest belongs to
         */
        virtual void interestRetransmit(const std::string &name,
                                     	unsigned int threadId) = 0;

        /**
         * Called by NDN-RTC every time Data segment has been received
         * @param name - 	Data segment name
         * @param threadId -ID of the media thread received Data segment belongs to
         * @param ndnPacketSize - 	Full size (including NDN packet overhead) of 
         *							Data segment packet (in  bytes)
         */
        virtual void dataReceived(const std::string &interestName,
                                  const std::string &dataName,
                                  unsigned int threadId,
                                  unsigned int ndnPacketSize) = 0;

        /**
         * EXPERIMENTAL API (11/26/2015)
         * !!! THIS CALL MUST BE USED INSTEAD OF dataReceived
         *
         * Called by NDN-RTC every time Data segment has been received
         * @param name          - Data segment name
         * @param threadId      - ID of the media thread received Data segment 
         *                      belongs to
         * @param payloadSize   - Payload size of the packet (segment size set 
         *                      by producer)
         * @param ndnPacketSize - Full size (including NDN packet overhead) of
         *						Data segment packet (in  bytes)
         * @param dataNonce     - nonce value of the Interest that retrieved this
         *                      data segment from consumer. Possible meanings:
         *                          1) value is 0: data was produced before any
         *                          Interest retrieved it
         *                          2) value equals interestNonce of previously 
         *                          issued Interest: data was retrieved from producer
         *                          and Dgen is valid
         *                          3) value is not 0 and doesn't equal interestNonce 
         *                          of any previous Interest: data was retrieved by 
         *                          another Consumer, dGen may be inaccurate
         * @param dGen          - generation delay for received data segment
         */
        virtual void dataReceivedX(const std::string &interestName,
                                   const std::string &dataName,
                                   unsigned int threadId,
                                   unsigned int payloadSize,
                                   unsigned int ndnPacketSize,
                                   int32_t dataNonce,
                                   int32_t dGen) = 0;
        
        /**
         * Called by NDN-RTC whenever any of the indicators has been updated
         * @param indicators - ARC indicators
         * @see ArcIndicators
         */
        virtual void updateIndicators(const ArcIndicators& indicators) = 0;

        /**
         * Called by NDN-RTC whenever any thread event occurred
         * @param event - Stream event
         * @see ThreadEvent
         */
        virtual void reportThreadEvent(const ThreadEvent& event) = 0;
    };
    
    /**
     * ARC's module callback (represents NDN-RTC)
     */
    class IRateAdaptationModuleCallback {
    	public:
    		/**
    		 * Called by ARC module whenever it wants to initiate challenge phase
    		 * @param threadId -	ID of the thread (usually, with higher bitrate) to 
    		 *						start fetching from
    		 * @param challengeLevel - 	Percentage (in the interval [0,1]) of thread 
    		 *							data to fetch
    		 */
    		virtual void onChallengePhaseStarted(unsigned int threadId, 
    											 double challengeLevel) = 0;
    		/**
    		 * Called by ARC module whenever it wants to stop Challenge phase
    		 */
    		virtual void onChallengePhaseStopped() = 0;

    		/**
    		 * Called by ARC whenever it wants to change challenging level.
    		 * NOTE: ARC module should initiate Challenge phase before calling 
    		 * this method
    		 * @param challengeLevel -	Percentate (in the interval [0,1]) of thread
    		 *							data to fetch
    		 * @see onChallengePhaseStarted
    		 */
    		virtual void onThreadChallenge(double challengeLevel) = 0;

    		/**
    		 * Called by ARC module when it decides that current network 
    		 * conditions are sufficient for fetching new thread (usually, 
    		 * with higher bandwidth)
             * @param threadId ID of the media thread wichi should be fetched
    		 * @see onChallengePhaseStarted
    		 */
			virtual void onThreadShouldSwitch(unsigned int threadId) = 0;
    };
}

#endif
