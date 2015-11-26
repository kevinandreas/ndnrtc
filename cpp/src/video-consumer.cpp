//
//  video-consumer.cpp
//  ndnrtc
//
//  Copyright 2013 Regents of the University of California
//  For licensing details see the LICENSE file.
//
//  Author:  Peter Gusev
//

#include "video-consumer.h"
#include "frame-buffer.h"
#include "pipeliner.h"
#include "chase-estimation.h"
#include "buffer-estimator.h"
#include "rtt-estimation.h"
#include "video-playout.h"
#include "interest-queue.h"
#include "video-decoder.h"
#include "arc-module.h"
#include "ndnrtc-namespace.h"

using namespace boost;
using namespace ndnlog;
using namespace ndnrtc;
using namespace ndnrtc::new_api;

//******************************************************************************
#pragma mark - construction/destruction
VideoConsumer::VideoConsumer(const GeneralParams& generalParams,
                             const GeneralConsumerParams& consumerParams,
                             IExternalRenderer* const externalRenderer):
Consumer(generalParams, consumerParams),
decoder_(new NdnVideoDecoder())
{
    setDescription("vconsumer");
    renderer_.reset(new ExternalVideoRendererAdaptor(externalRenderer));
    decoder_->setFrameConsumer(getRenderer().get());
}

VideoConsumer::~VideoConsumer()
{
    
}

//******************************************************************************
#pragma mark - public
int
VideoConsumer::init(const ConsumerSettings& settings,
                    const std::string& threadName)
{
    int res = RESULT_OK;
    
    LogInfoC << "unix timestamp: " << std::fixed << std::setprecision(6)
    << NdnRtcUtils::unixTimestamp() << std::endl;

    if (RESULT_GOOD(Consumer::init(settings, threadName)))
    {
        initArc();
        interestQueue_->registerCallback(this);
        
        pipeliner_->setUseKeyNamespace(true);
        pipeliner_->initialize();
        
        decoder_->init(((VideoThreadParams*)getCurrentThreadParameters())->coderParams_);
        
        playout_.reset(new VideoPlayout(this, statStorage_));
        playout_->setLogger(logger_);
        playout_->init(decoder_.get());
        ((VideoPlayout*)playout_.get())->onFrameSkipped_ = boost::bind(&VideoConsumer::onFrameSkipped, this, _1, _2, _3, _4, _5);
        
        LogInfoC << "initialized" << std::endl;
        
        return res;
    }
    
    return RESULT_ERR;
}

int
VideoConsumer::start()
{
    if (RESULT_GOOD(Consumer::start()))
        LogInfoC << "started" << std::endl;
    else
        return RESULT_ERR;
    
    return RESULT_OK;
}

int
VideoConsumer::stop()
{
    if (RESULT_GOOD(Consumer::stop()))
        LogInfoC << "stopped" << std::endl;
    else
        return RESULT_ERR;
    
    return RESULT_OK;
}

void
VideoConsumer::setLogger(ndnlog::new_api::Logger *logger)
{
    getRenderer()->setLogger(logger);
    decoder_->setLogger(logger);
    
    Consumer::setLogger(logger);
}

void
VideoConsumer::onStateChanged(const int& oldState, const int& newState)
{
    Consumer::onStateChanged(oldState, newState);
#if 0
    if (rateControl_.get())
    {
        if (newState == Pipeliner::StateFetching)
            rateControl_->start();
        
        if (oldState == Pipeliner::StateFetching)
            rateControl_->stop();
    }
#endif
}

void
VideoConsumer::playbackEventOccurred(PlaybackEvent event,
                                     unsigned int frameSeqNo)
{
    if (observer_)
    {
        lock_guard<mutex> scopedLock(observerMutex_);
        observer_->onPlaybackEventOccurred(event, frameSeqNo);
    }
}

void
VideoConsumer::triggerRebuffering()
{
    Consumer::triggerRebuffering();
    decoder_->reset();
}

//******************************************************************************
void
VideoConsumer::onFrameSkipped(PacketNumber playbackNo, PacketNumber sequenceNo,
                              PacketNumber pairedNo, bool isKey,
                              double assembledLevel)
{
    // empty
}

void
VideoConsumer::onTimeout(const shared_ptr<const Interest>& interest)
{
#if 0
    if (rateControl_.get())
        rateControl_->interestTimeout(interest);
#endif
}

void
VideoConsumer::initArc()
{
    for (int i = 0; i < settings_.streamParams_.mediaThreads_.size(); i++)
    {
        VideoThreadParams *threadParams = (VideoThreadParams*)settings_.streamParams_.mediaThreads_[i];
        ArcModule::ThreadEntry entry = {i, (double)threadParams->coderParams_.startBitrate_, 0.2};
        arcThreadArray_.push_back(entry);
    }
     
    arcModule_.reset(new ArcModule(NdnRtcUtils::getIoService()));
    arcModule_->initialize(this, IRateAdaptationModule::CodecMode::CodecModeNormal, arcThreadArray_);
}

void
VideoConsumer::onChallengePhaseStarted(unsigned int threadId,
                                       double challengeLevel)
{
    
    std::string threadName = settings_.streamParams_.mediaThreads_[threadId]->threadName_;
    FrameSegmentsInfo segInfo = settings_.streamParams_.mediaThreads_[threadId]->getSegmentsInfo();
    
    pipeliner_->startChallenging(threadName, challengeLevel, segInfo);
}

void
VideoConsumer::onChallengePhaseStopped()
{
    pipeliner_->stopChallenging();
}

void
VideoConsumer::onThreadChallenge(double challengeLevel)
{
    pipeliner_->newChallengeLevel(challengeLevel);
}

void
VideoConsumer::onThreadShouldSwitch(unsigned int threadId)
{
    currentThreadIdx_ = threadId;
    
    LogInfoC << "ARC: switching to " << getCurrentThreadName()
    << " initiated" << std::endl;
    
    pipeliner_->threadSwitched();
    
    if (observer_)
    {
        lock_guard<mutex> scopedLock(observerMutex_);
        observer_->onThreadSwitched(getCurrentThreadName());
    }
}

void
VideoConsumer::onInterestIssued(const boost::shared_ptr<const ndn::Interest>& interest)
{
    int threadIdx = threadIdxMap_[NdnRtcNamespace::getThreadName(interest->getName())];
    uint32_t nonce = NdnRtcUtils::blobToNonce(interest->getNonce());
    
    arcModule_->interestExpressed(interest->getName().toUri(),threadIdx, nonce);
}

void
VideoConsumer::onDataArrived(const boost::shared_ptr<const Interest>& interest,
                             const boost::shared_ptr<ndn::Data>& data)
{
    std::string threadName = NdnRtcNamespace::getThreadName(data->getName());
    SegmentData segmentData;
    
    if (RESULT_GOOD(SegmentData::segmentDataFromRaw(data->getContent().size(),
                                                    data->getContent().buf(),
                                                    segmentData)))
    {
        SegmentData::SegmentMetaInfo segmentMeta = segmentData.getMetadata();
        
        arcModule_->dataReceivedX(interest->getName().toUri(),
                                 data->getName().toUri(),
                                 getThreadIdx(threadName),
                                 data->getDefaultWireEncoding().size(),
                                 segmentMeta.interestNonce_,
                                 segmentMeta.generationDelayMs_);
    }
}