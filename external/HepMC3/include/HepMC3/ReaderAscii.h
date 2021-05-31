// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_READERASCII_H
#define HEPMC3_READERASCII_H
///
/// @file  ReaderAscii.h
/// @brief Definition of class \b ReaderAscii
///
/// @class HepMC3::ReaderAscii
/// @brief GenEvent I/O parsing for structured text files
///
/// @ingroup IO
///
#include <set>
#include <string>
#include <fstream>
#include <istream>
#include <iterator>
#include "HepMC3/Reader.h"
#include "HepMC3/GenEvent.h"


namespace HepMC3 {


class ReaderAscii : public Reader {
public:

    /// @brief Constructor
    ReaderAscii(const std::string& filename);
    /// The ctor to read from stdin
    ReaderAscii(std::istream &);
    /// @brief Destructor
    ~ReaderAscii();

    /// @brief skip events
    bool skip(const int)  override;
    
    /// @brief Load event from file
    ///
    /// @param[out] evt Event to be filled
    bool read_event(GenEvent& evt)  override;

    /// @todo No-arg version returning GenEvent?

    /// @brief Return status of the stream
    bool failed()  override;

    /// @todo Implicit cast to bool = !failed()?

    /// @brief Close file stream
    void close()  override;

private:

    /// @brief Unsecape '\' and '\n' characters in string
    std::string unescape(const std::string& s);

    /// @name Read helpers
    //@{

    /// @brief Parse event
    ///
    /// Helper routine for parsing event information
    /// @param[out] evt Event that will be filled with new data
    /// @param[in]  buf Line of text that needs to be parsed
    /// @return vertices count and particles count for verification
    std::pair<int,int> parse_event_information(GenEvent &evt, const char *buf);

    /// @brief Parse weight value lines
    ///
    /// Helper routine for parsing weight value information
    /// @param[out] evt Event whose GenWeights will be filled with weight values
    /// @param[in]  buf Line of text that needs to be parsed
    ///
    bool parse_weight_values(GenEvent &evt, const char *buf);

    /// @brief Parse units
    ///
    /// Helper routine for parsing units information
    /// @param[out] evt Event that will be filled with unit information
    /// @param[in]  buf Line of text that needs to be parsed
    ///
    bool parse_units(GenEvent &evt, const char *buf);

    /// @brief Parse struct GenPdfInfo information
    ///
    /// Helper routine for parsing PDF information
    /// @param[out] evt Event that will be filled with unit information
    /// @param[in]  buf Line of text that needs to be parsed
    bool parse_pdf_info(GenEvent &evt, const char *buf);

    /// @brief Parse struct GenHeavyIon information
    ///
    /// Helper routine for parsing heavy ion information
    /// @param[out] evt Event that will be filled with unit information
    /// @param[in]  buf Line of text that needs to be parsed
    bool parse_heavy_ion(GenEvent &evt, const char *buf);

    /// @brief Parse struct GenCrossSection information
    ///
    /// Helper routine for parsing cross-section information
    /// @param[out] evt Event that will be filled with unit information
    /// @param[in]  buf Line of text that needs to be parsed
    bool parse_cross_section(GenEvent &evt, const char *buf);

    /// @brief Parse vertex
    ///
    /// Helper routine for parsing single event information
    /// @param[out] evt Event that will contain parsed vertex
    /// @param[in] buf Line of text that needs to be parsed
    ///
    bool parse_vertex_information(GenEvent &evt, const char *buf);

    /// @brief Parse particle
    ///
    /// Helper routine for parsing single particle information
    /// @param[out] evt Event that will contain parsed particle
    /// @param[in] buf Line of text that needs to be parsed
    bool parse_particle_information(GenEvent &evt, const char *buf);

    /// @brief Parse attribute
    ///
    /// Helper routine for parsing single attribute information
    /// @param[out] evt Event that will contain parsed attribute
    /// @param[in] buf Line of text that needs to be parsed
    bool parse_attribute(GenEvent &evt, const char *buf);

    /// @brief Parse run-level attribute.
    ///
    /// Helper routine for parsing single attribute information
    /// @param[in] buf Line of text that needs to be parsed
    bool parse_run_attribute(const char *buf);

    /// @brief Parse run-level weight names.
    ///
    /// Helper routine for parsing a line with information about
    /// weight names.
    /// @param[in] buf Line of text that needs to be parsed
    bool parse_weight_names(const char *buf);

    /// @brief Parse run-level tool information.
    ///
    /// Helper routine for parsing a line with information about
    /// tools being used.
    /// @param[in] buf Line of text that needs to be parsed
    bool parse_tool(const char *buf);
    //@}


private:

    std::ifstream m_file; //!< Input file
    std::istream* m_stream; ///< For ctor when reading from stdin
    bool m_isstream; ///< toggles usage of m_file or m_stream


    /** @brief Store attributes global to the run being written/read. */
    std::map< std::string, std::shared_ptr<Attribute> > m_global_attributes;

    /** @brief Temp storage for  outgoing particle ids */
    std::map<GenVertexPtr, std::set<int> >  m_forward_mothers;
    /** @brief Temp storage for  prod vertex ids */
    std::map<GenParticlePtr, int >  m_forward_daughters;

};


} // namespace HepMC3

#endif
