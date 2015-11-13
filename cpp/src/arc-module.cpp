//
//  arc-module.cpp
//  ndnrtc
//
//  Created by Takahiro Yoneda on 11/13/15.
//  Copyright 2013-2015 Regents of the University of California
//

#include "arc-module.h"

using namespace ndnrtc;

int ArcModule::initialize(IRateAdaptationModuleCallback* const callback,
						const CodecMode& codecMode,
        				std::vector<ThreadEntry> mediaThreads)
{
	callback_ = callback;
	return 0;
}
        
void interestExpressed(const std::string &name,
                               unsigned int threadId)
{}

void interestRetransmit(const std::string &name,
                             	unsigned int threadId)
{}

void dataReceived(const std::string &name,
                          unsigned int threadId,
                          unsigned int ndnPacketSize)
{}

void updateIndicators(const ArcModule::ArcIndicators& indicators)
{}

void reportThreadEvent(const ArcModule::ThreadEvent& event)
{}