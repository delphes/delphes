// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_READERROOTTREE_H
#define HEPMC3_READERROOTTREE_H
/**
 *  @file  ReaderRootTree.h
 *  @brief Definition of \b class ReaderRootTree
 *
 *  @class HepMC3::ReaderRootTree
 *  @brief GenEvent I/O parsing and serialization for root files  based on root TTree
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
#include "TTree.h"
#include "TBranch.h"

namespace HepMC3
{

class ReaderRootTree : public Reader
{
//
// Constructors
//
public:
    /** @brief Default constructor */
    ReaderRootTree(const std::string &filename);
    /** @brief Constructor with tree name*/
    ReaderRootTree(const std::string &filename, const std::string &treename, const std::string &branchname);

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
    bool read_event(GenEvent &evt)   override;

    /** @brief Close file */
    void close()  override;

    /** @brief Get file  error state */
    bool failed()  override;

private:
    /** @brief init routine */
    bool init();
//
// Fields
//
private:
    TFile* m_file;         //!< File handler
public:
    TTree* m_tree;//!< Tree handler. Public to allow simple access, e.g. custom branches.
private:
    int   m_events_count; //!< Events count. Needed to read the tree
    GenEventData* m_event_data; //!< Pointer to structure that holds event data
    GenRunInfoData* m_run_info_data; //!< Pointer to structure that holds run info data
    std::string m_tree_name; //!< Name of TTree
    std::string m_branch_name; //!< Name of TBranch in TTree
};

} // namespace HepMC3

#endif
