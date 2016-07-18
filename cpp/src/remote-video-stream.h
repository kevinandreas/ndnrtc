// 
// remote-video-stream.h
//
//  Created by Peter Gusev on 29 June 2016.
//  Copyright 2013-2016 Regents of the University of California
//

#ifndef __remote_video_stream_h__
#define __remote_video_stream_h__

#include "remote-stream-impl.h"

namespace ndnrtc{
	class IPlayoutControl;
	class VideoPlayout;
	class PipelineControl;

	class RemoteVideoStreamImpl : public RemoteStreamImpl
	{
	public:
		RemoteVideoStreamImpl(boost::asio::io_service& io, 
			const boost::shared_ptr<ndn::Face>& face,
			const boost::shared_ptr<ndn::KeyChain>& keyChain,
			const std::string& streamPrefix);
		
		void initiateFetching();

	private:
	};
}

#endif