// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_WRITERROOTTREE_H
#define HEPMC3_WRITERROOTTREE_H
/**
 *  @file  WriterRootTree.h
 *  @brief Definition of \b class WriterRootTree
 *
 *  @class HepMC3::WriterRootTree
 *  @brief GenEvent I/O serialization for root files based on root TTree
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
#include "TTree.h"
#else
class TFile;
class TTree;
#endif

namespace HepMC3
{
class WriterRootTree : public Writer
{
//
// Constructors
//
public:
    /** @brief Default constructor
     *  @warning If file exists, it will be overwritten
     */
    WriterRootTree(const std::string &filename,
                   std::shared_ptr<GenRunInfo> run = std::shared_ptr<GenRunInfo>());
    /** @brief Constructor with tree name*/
    WriterRootTree(const std::string &filename, const std::string &treename, const std::string &branchname,
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
    bool failed()  override;

private:
    /** @brief init routine */
    bool init(std::shared_ptr<GenRunInfo> run);
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
    std::string m_tree_name;//!< Name of TTree
    std::string m_branch_name; //!< Name of TBranch in TTree
};

} // namespace HepMC3

#endif
