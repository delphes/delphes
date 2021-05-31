// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file WriterRootTree.cc
 *  @brief Implementation of \b class WriterRootTree
 *
 */
#include <cstdio>  // sprintf
#include "HepMC3/WriterRootTree.h"
#include "HepMC3/Version.h"
// ROOT header files
#include "TFile.h"
#include "TTree.h"

namespace HepMC3
{
HEPMC3_DECLARE_WRITER_FILE(WriterRootTree)

WriterRootTree::WriterRootTree(const std::string &filename, std::shared_ptr<GenRunInfo> run):
    m_tree(0),
    m_events_count(0),
    m_tree_name("hepmc3_tree"),
    m_branch_name("hepmc3_event")
{
    m_file = TFile::Open(filename.c_str(), "RECREATE");
    if (!init(run)) return;
}

WriterRootTree::WriterRootTree(const std::string &filename, const std::string &treename, const std::string &branchname, std::shared_ptr<GenRunInfo> run):
    m_tree(0),
    m_events_count(0),
    m_tree_name(treename.c_str()),
    m_branch_name(branchname.c_str())
{
    m_file = TFile::Open(filename.c_str(), "RECREATE");
    if (!init(run)) return;
}

bool WriterRootTree::init(std::shared_ptr<GenRunInfo> run )
{
    if ( !m_file->IsOpen() )
    {
        HEPMC3_ERROR("WriterRootTree: problem opening file: " << m_file->GetName())
        return false;
    }
    m_event_data = new GenEventData();
    m_run_info_data = new GenRunInfoData();
    set_run_info(run);
    if ( run_info() ) run_info()->write_data(*m_run_info_data);
    m_tree = new TTree(m_tree_name.c_str(), "hepmc3_tree");
    m_tree->Branch(m_branch_name.c_str(), m_event_data);
    m_tree->Branch("GenRunInfo", m_run_info_data);
    return true;
}

void WriterRootTree::write_event(const GenEvent &evt)
{
    if ( !m_file->IsOpen() ) return;
    bool refill = false;
    if ( evt.run_info()&&(!run_info() || (run_info() != evt.run_info())))  { set_run_info(evt.run_info()); refill = true;}
    if (refill)
    {
        m_run_info_data->weight_names.clear();
        m_run_info_data->tool_name.clear();
        m_run_info_data->tool_version.clear();
        m_run_info_data->tool_description.clear();
        m_run_info_data->attribute_name.clear();
        m_run_info_data->attribute_string.clear();
        run_info()->write_data(*m_run_info_data);
    }



    m_event_data->particles.clear();
    m_event_data->vertices.clear();
    m_event_data->links1.clear();
    m_event_data->links2.clear();
    m_event_data->attribute_id.clear();
    m_event_data->attribute_name.clear();
    m_event_data->attribute_string.clear();

    evt.write_data(*m_event_data);
    m_tree->Fill();
    ++m_events_count;
}


void WriterRootTree::write_run_info() {}

void WriterRootTree::close()
{
    m_file->WriteTObject(m_tree);
    m_file->Close();
    delete m_event_data;
    delete m_run_info_data;
}

bool WriterRootTree::failed()
{
    if ( !m_file->IsOpen() ) return true;

    return false;
}

} // namespace HepMC3
