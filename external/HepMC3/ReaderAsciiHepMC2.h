// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_READER_ASCII_HEPMC2_H
#define HEPMC3_READER_ASCII_HEPMC2_H
/**
 *  @file  ReaderAsciiHepMC2.h
 *  @brief Definition of \b class ReaderAsciiHepMC2
 *
 *  @class HepMC3::ReaderAsciiHepMC2
 *  @brief Parser for HepMC2 I/O files
 *
 *  @ingroup IO
 *
 */
#include "HepMC3/Reader.h"

#include "HepMC3/GenEvent.h"

#include <string>
#include <fstream>
#include <istream>

namespace HepMC3 {



class ReaderAsciiHepMC2 : public Reader {
//
// Constructors
//
public:
    /** @brief Default constructor */
    ReaderAsciiHepMC2(const std::string& filename);

    /// The ctor to read from stdin
    ReaderAsciiHepMC2(std::istream &);

    /// @brief Destructor
    ~ReaderAsciiHepMC2();
//
// Functions
//
public:
    /// @brief skip events
    bool skip(const int)  override;

    /** @brief Implementation of Reader::read_event */
    bool read_event(GenEvent &evt)  override;

    /// @brief Return status of the stream
    bool failed()  override;

    /// @brief Close file stream
    void close()  override;

private:
    /** @brief Parse event
     *
     *  Helper routine for parsing event information
     *  @param[out] evt Event that will be filled with new data
     *  @param[in]  buf Line of text that needs to be parsed
     */
    int parse_event_information(GenEvent &evt, const char *buf);

    /** @brief Parse units
     *
     *  Helper routine for parsing unit information
     *  @param[out] evt Event that will be filled with unit information
     *  @param[in]  buf Line of text that needs to be parsed
     */
    bool parse_units(GenEvent &evt, const char *buf);

    /** @brief Parse vertex
     *
     *  Helper routine for parsing single event information
     *  @param[in] buf Line of text that needs to be parsed
     */
    int parse_vertex_information(const char *buf);

    /** @brief Parse particle
     *
     *  Helper routine for parsing single particle information
     *  @param[in] buf Line of text that needs to be parsed
     */
    int parse_particle_information(const char *buf);

    /** @brief Parse weight names
     *
     *  Helper routine for parsing weight names
     *  @param[in] buf Line of text that needs to be parsed
     */
    bool parse_weight_names(const char *buf);

    /** @brief Parse heavy ion information
     *
     *  Helper routine for parsing heavy ion information
     *  @param[out] evt Event that will be filled with new data
     *  @param[in]  buf Line of text that needs to be parsed
     */
    bool parse_heavy_ion(GenEvent &evt, const char *buf);

    /** @brief Parse pdf information
     *
     *  Helper routine for parsing pdf information
     *  @param[out] evt Event that will be filled with new data
     *  @param[in]  buf Line of text that needs to be parsed
     */
    bool parse_pdf_info(GenEvent &evt, const char *buf);


    /** @brief Parse pdf information
    *
    *  Helper routine for parsing cross-section information
    *  @param[out] evt Event that will be filled with new data
    *  @param[in]  buf Line of text that needs to be parsed
    */
    bool parse_xs_info(GenEvent &evt, const char *buf);



//
// Fields
//
private:
    std::ifstream m_file; //!< Input file
    std::istream* m_stream; ///< For ctor when reading from stdin
    bool m_isstream; ///< toggles usage of m_file or m_stream

    std::vector<GenVertexPtr>   m_vertex_cache;        //!< Vertex cache
    std::vector<int>            m_vertex_barcodes;     //!< Old vertex barcodes

    std::vector<GenParticlePtr> m_particle_cache;      //!< Particle cache
    std::vector<int>            m_end_vertex_barcodes; //!< Old end vertex barcodes

    GenEvent*              m_event_ghost;                      //!< To save particle and verstex attributes.
    std::vector<GenParticlePtr> m_particle_cache_ghost;//!< Particle cache for attributes
    std::vector<GenVertexPtr>   m_vertex_cache_ghost;        //!< Vertex cache for attributes
};

} // namespace HepMC3

#endif
