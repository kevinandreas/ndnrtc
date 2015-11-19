//
//  arc-module.cpp
//  ndnrtc
//
//  Created by Takahiro Yoneda on 11/13/15.
//  Copyright 2013-2015 Regents of the University of California
//

#include "arc-module.h"

using namespace ndnrtc;

const ArcModule::ArcIndicators ArcModule::ZeroIndicators = {0,0.,0.,0.,0.,0.,0.,0.,ConsumerPhaseInactive};

int ArcModule::initialize(IRateAdaptationModuleCallback* const callback,
						const CodecMode& codecMode,
        				std::vector<ThreadEntry> mediaThreads)
{
	callback_ = callback;
	return 0;
}
        
void ArcModule::interestExpressed(const std::string &name,
                               unsigned int threadId)
{}

void ArcModule::interestRetransmit(const std::string &name,
                             	unsigned int threadId)
{}

void ArcModule::dataReceived(const std::string &name,
                          unsigned int threadId,
                          unsigned int ndnPacketSize)
{}

void ArcModule::updateIndicators(const ArcModule::ArcIndicators& indicators)
{}

void ArcModule::reportThreadEvent(const ArcModule::ThreadEvent& event)
{}