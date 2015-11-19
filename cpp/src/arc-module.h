//
//  arc-module.h
//  ndnrtc
//
//  Created by Takahiro Yoneda on 11/13/15.
//  Copyright 2013-2015 Regents of the University of California
//

#include "arc-interface.h"

namespace ndnrtc {
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
        
        ArcModule(boost::asio::io_service& ioService):ArcModuleIoBase(ioService){}
        ~ArcModule(){}
        
        int initialize(IRateAdaptationModuleCallback* const callback,
                       const CodecMode& codecMode,
                       std::vector<ThreadEntry> mediaThreads);
        void interestExpressed(const std::string &name,
                               unsigned int threadId);
        void interestRetransmit(const std::string &name,
                                unsigned int threadId);
        void dataReceived(const std::string &name,
                          unsigned int threadId,
                          unsigned int ndnPacketSize);
        void updateIndicators(const ArcIndicators& indicators);
        void reportThreadEvent(const ThreadEvent& event);
        
    private:
        IRateAdaptationModuleCallback* callback_;
    };
};