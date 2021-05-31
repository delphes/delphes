// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_WRITERROOT_H
#define HEPMC3_WRITERROOT_H
/**
 *  @file  WriterRoot.h
 *  @brief Definition of \b class WriterRoot
 *
 *  @class HepMC3::WriterRoot
 *  @brief GenEvent I/O serialization for root files
 *
 *  If HepMC was compiled with path to ROOT available, this class can be used
 *  for root writing in the same manner as with HepMC::WriterAscii class.
 *
 *  @ingroup IO
 *
 */
#include <string>
#include <memory>
#include "HepMC3/Writer.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/Data/GenEventData.h"
#include "HepMC3/Data/GenRunInfoData.h"

// ROOT header files
#ifdef __CINT__
#include "TFile.h"
#else
class TFile;
#endif

namespace HepMC3 {

class WriterRoot : public Writer {
//
// Constructors
//
public:
    /** @brief Default constructor
     *  @warning If file exists, it will be overwritten
     */
    WriterRoot(const std::string& filename,
               std::shared_ptr<GenRunInfo> run = std::shared_ptr<GenRunInfo>());

//
// Functions
//
public:
    /** @brief Write event to file
     *
     *  @param[in] evt Event to be serialized
     */
    void write_event(const GenEvent &evt) override;

    /** @brief Write the GenRunInfo object to file. */
    void write_run_info();

    /** @brief Close file stream */
    void close() override;

    /** @brief Get stream error state flag */
    bool failed() override;
//
// Fields
//
private:
    TFile* m_file;         //!< File handler
    int    m_events_count; //!< Events count. Needed to generate unique object name
};

} // namespace HepMC3

#endif
