// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_READERROOT_H
#define HEPMC3_READERROOT_H
/**
 *  @file  ReaderRoot.h
 *  @brief Definition of \b class ReaderRoot
 *
 *  @class HepMC3::ReaderRoot
 *  @brief GenEvent I/O parsing and serialization for root files
 *
 *  If HepMC was compiled with path to ROOT available, this class can be used
 *  for root file I/O in the same manner as with HepMC::ReaderAscii class.
 *
 *  @ingroup IO
 *
 */
#include <string>
#include "HepMC3/Reader.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/Data/GenEventData.h"
#include "HepMC3/Data/GenRunInfoData.h"

// ROOT header files
#include "TFile.h"
#include "TKey.h"

namespace HepMC3 {

class ReaderRoot : public Reader {
//
// Constructors
//
public:
    /** @brief Default constructor */
    ReaderRoot(const std::string &filename);

//
// Functions
//
public:
    /// @brief skip events
    bool skip(const int)  override;

    /** @brief Read event from file
     *
     *  @param[out] evt Contains parsed event
     */
    bool read_event(GenEvent &evt)  override;

    /** @brief Close file stream */
    void close()  override;

    /** @brief Get stream error state */
    bool failed()  override;
//
// Fields
//
private:
    TFile* m_file; //!< File handler
    TIter* m_next; //!< Iterator for event reading
};

} // namespace HepMC3

#endif
