#ifndef __FASTJET_CONFIG_H__
#define __FASTJET_CONFIG_H__

// by default, use an automatically generated config_auto.h
// unless it's a windows machine in which case
#ifndef WIN32
#include "config_auto.h"
#else 
#include "fastjet/config_win.h"
#endif // WIN32

#endif // __FASTJET_CONFIG_H__
