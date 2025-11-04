#ifndef __FASTJET_CONFIG_H__
#define __FASTJET_CONFIG_H__

// Required definitions for exporting static variables in windows builds (for DLLs)
// This is only needed for static data variables since we use
// the CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON automation when building with cmake.
// That automation handles all member functions.
// So, when making a static variable please add in the beginning of a variable
// definition, like a keyword. It is very important to use the visibility relevant
// to the library you are working in, see below for possibilities!
// e.g.: FASTJET_WINDLL static bool verbosity; 
// inside a SomeClass.h, for instance
// Similarly for declarations you must prefix the appropriate WINDLL
// e.g. FASTJET_WINDLL static bool SomeClass::verbosity = true; // in SomeClass.(hh|cc)
#ifdef _WIN32
    #if defined(fastjet_EXPORTS)
        #define FASTJET_WINDLL __declspec(dllexport) // Export when building the DLL
    #else
        #define FASTJET_WINDLL __declspec(dllimport) // Import when using the DLL
    #endif

    #if defined(fastjettools_EXPORTS)
        #define FASTJET_TOOLS_WINDLL __declspec(dllexport) // Export when building the DLL
    #else
        #define FASTJET_TOOLS_WINDLL __declspec(dllimport) // Import when using the DLL
    #endif

    #if defined(fastjetplugins_EXPORTS)
        #define FASTJET_PLUGINS_WINDLL __declspec(dllexport) // Export when building the DLL
    #else
        #define FASTJET_PLUGINS_WINDLL __declspec(dllimport) // Import when using the DLL
    #endif
#else
    // For Linux/macOS
    #define FASTJET_WINDLL
    #define FASTJET_TOOLS_WINDLL
    #define FASTJET_PLUGINS_WINDLL
#endif

// by default, use an automatically generated config_auto.h
// unless it's a windows machine in which case
#include "fastjet/config_auto.h"

#endif // __FASTJET_CONFIG_H__
