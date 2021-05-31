// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_READER_H
#define HEPMC3_READER_H
///
/// @file  Reader.h
/// @brief Definition of interface \b Reader
///
/// @class HepMC3::Reader
/// @brief Base class for all I/O readers
///
/// @ingroup IO
///

#include "HepMC3/GenRunInfo.h"

namespace HepMC3 {

// Forward declaration
class GenEvent;

class Reader {
public:
    ///Constructor
    Reader() {}

    /// Virtual destructor
    virtual ~Reader() {}

    /// skip or fast forward reading of some events
    virtual bool skip(const int) { return !failed();}

    /// Fill next event from input into @a evt
    virtual bool read_event(GenEvent& evt) = 0;
    /** @brief Get file and/or stream error state */
    virtual bool failed()=0;
    /** @brief Close file and/or stream */
    virtual void close()=0;

    /// Get the global GenRunInfo object.
    std::shared_ptr<GenRunInfo> run_info() const {
        return m_run_info;
    }

///deleted copy constructor
    Reader(const Reader&) = delete; 
///deleted copy assignment operator
    Reader& operator = (const Reader &) = delete;            
    /// Set options
    void set_options(const std::map<std::string, std::string>& options)
    {
    m_options=options;
    }
    /// Set options
    std::map<std::string, std::string> get_options() const
    {
    return m_options;
    }
protected:
    /// Set the global GenRunInfo object.
    void set_run_info(std::shared_ptr<GenRunInfo> run) {
        m_run_info = run;
    }
        /// options
    std::map<std::string, std::string> m_options;
private:
    /// The global GenRunInfo object.
    std::shared_ptr<GenRunInfo> m_run_info;
};


} // namespace HepMC3

#endif
