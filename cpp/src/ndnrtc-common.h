//
//  ndnrtc_common.h
//  ndnrtc
//
//  Copyright 2013 Regents of the University of California
//  For licensing details see the LICENSE file.
//
//  Author:  Peter Gusev 
//  Created: 8/7/13
//

#ifndef ndnrtc_ndnrtc_common_h
#define ndnrtc_ndnrtc_common_h

//#define DEBUG
//#define NDN_LOGGING
//#define NDN_DETAILED
#define NDN_TRACE
#define NDN_INFO
#define NDN_WARN
#define NDN_ERROR
#define NDN_DEBUG

#include <stdexcept>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <cstring>

#include "ndnlib.h"
#include "webrtc.h"
#include "simple-log.h"

using namespace std;
using namespace ndn;
using namespace ptr_lib;
using namespace ndnlog;


#endif