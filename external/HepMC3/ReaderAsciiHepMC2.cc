// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file ReaderAsciiHepMC2.cc
 *  @brief Implementation of \b class ReaderAsciiHepMC2
 *
 */
#include "HepMC3/ReaderAsciiHepMC2.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenHeavyIon.h"
#include "HepMC3/GenPdfInfo.h"
#include "HepMC3/Setup.h"
#include <cstring>
#include <cstdlib>
namespace HepMC3 {

ReaderAsciiHepMC2::ReaderAsciiHepMC2(const std::string& filename):
    m_file(filename), m_stream(0), m_isstream(false) {
    if( !m_file.is_open() ) {
        HEPMC3_ERROR( "ReaderAsciiHepMC2: could not open input file: "<<filename )
    }
    set_run_info(std::make_shared<GenRunInfo>());
    m_event_ghost= new GenEvent();
}
// Ctor for reading from stdin
ReaderAsciiHepMC2::ReaderAsciiHepMC2(std::istream & stream)
    : m_stream(&stream), m_isstream(true)
{
    if( !m_stream->good() ) {
        HEPMC3_ERROR( "ReaderAsciiHepMC2: could not open input stream " )
    }
    set_run_info(std::make_shared<GenRunInfo>());
    m_event_ghost= new GenEvent();
}

ReaderAsciiHepMC2::~ReaderAsciiHepMC2() { if (m_event_ghost) { m_event_ghost->clear(); delete m_event_ghost; m_event_ghost=nullptr; } if (!m_isstream) close(); }

bool ReaderAsciiHepMC2::skip(const int n)
{
    const size_t       max_buffer_size=512*512;
    char               buf[max_buffer_size];
    int nn=n;
    while(!failed()) {
        char peek;
        if ( (!m_file.is_open()) && (!m_isstream) ) return false;
        m_isstream ? peek = m_stream->peek() : peek = m_file.peek();
        if( peek=='E' ) nn--;
        if (nn<0) return true;
        m_isstream ? m_stream->getline(buf,max_buffer_size) : m_file.getline(buf,max_buffer_size);
    }
    return true;
}

bool ReaderAsciiHepMC2::read_event(GenEvent &evt) {
    if ( (!m_file.is_open()) && (!m_isstream) ) return false;

    char               peek;
    const size_t  max_buffer_size=512*512;
    const size_t  max_weights_size=256;
    char          buf[max_buffer_size];
    bool          parsed_event_header            = false;
    bool          is_parsing_successful          = true;
    int           parsing_result                 = 0;
    unsigned int  vertices_count                 = 0;
    unsigned int  current_vertex_particles_count = 0;
    unsigned int  current_vertex_particles_parsed= 0;

    evt.clear();
    evt.set_run_info(run_info());

    // Empty cache
    m_vertex_cache.clear();
    m_vertex_barcodes.clear();

    m_particle_cache.clear();
    m_end_vertex_barcodes.clear();
    m_particle_cache_ghost.clear();
    //
    // Parse event, vertex and particle information
    //
    while(!failed()) {
        m_isstream ? m_stream->getline(buf,max_buffer_size) : m_file.getline(buf,max_buffer_size);
        if( strlen(buf) == 0 ) continue;
        // Check for IO_GenEvent header/footer
        if( strncmp(buf,"HepMC",5) == 0 ) {
            if( strncmp(buf,"HepMC::Version",14) != 0 && strncmp(buf,"HepMC::IO_GenEvent",18)!=0 )
            {
                HEPMC3_WARNING( "ReaderAsciiHepMC2: found unsupported expression in header. Will close the input." )
                std::cout<<buf<<std::endl;
                m_isstream ? m_stream->clear(std::ios::eofbit) : m_file.clear(std::ios::eofbit);
            }
            if(parsed_event_header) {
                is_parsing_successful = true;
                break;
            }
            continue;
        }
        switch(buf[0]) {
        case 'E':
            parsing_result = parse_event_information(evt,buf);
            if(parsing_result<0) {
                is_parsing_successful = false;
                HEPMC3_ERROR( "ReaderAsciiHepMC2: HEPMC3_ERROR parsing event information" )
            }
            else {
                vertices_count = parsing_result;
                m_vertex_cache.reserve(vertices_count);
                m_particle_cache.reserve(vertices_count*3);
                m_vertex_barcodes.reserve(vertices_count);
                m_end_vertex_barcodes.reserve(vertices_count*3);
                is_parsing_successful = true;
            }
            parsed_event_header = true;
            break;
        case 'V':
            // If starting new vertex: verify if previous was fully parsed

            /** @bug HepMC2 files produced with Pythia8 are known to have wrong
                     information about number of particles in vertex. Hence '<' sign */
            if(current_vertex_particles_parsed < current_vertex_particles_count) {
                is_parsing_successful = false;
                break;
            }
            current_vertex_particles_parsed = 0;

            parsing_result = parse_vertex_information(buf);

            if(parsing_result<0) {
                is_parsing_successful = false;
                HEPMC3_ERROR( "ReaderAsciiHepMC2: HEPMC3_ERROR parsing vertex information" )
            }
            else {
                current_vertex_particles_count = parsing_result;
                is_parsing_successful = true;
            }
            break;
        case 'P':

            parsing_result   = parse_particle_information(buf);

            if(parsing_result<0) {
                is_parsing_successful = false;
                HEPMC3_ERROR( "ReaderAsciiHepMC2: HEPMC3_ERROR parsing particle information" )
            }
            else {
                ++current_vertex_particles_parsed;
                is_parsing_successful = true;
            }
            break;
        case 'U':
            is_parsing_successful = parse_units(evt,buf);
            break;
        case 'F':
            is_parsing_successful = parse_pdf_info(evt,buf);
            break;
        case 'H':
            is_parsing_successful = parse_heavy_ion(evt,buf);
            break;
        case 'N':
            is_parsing_successful = parse_weight_names(buf);
            break;
        case 'C':
            is_parsing_successful = parse_xs_info(evt,buf);
            break;
        default:
            HEPMC3_WARNING( "ReaderAsciiHepMC2: skipping unrecognised prefix: " << buf[0] )
            is_parsing_successful = true;
            break;
        }

        if( !is_parsing_successful ) break;

        // Check for next event
        m_isstream ? peek = m_stream->peek() : peek = m_file.peek();
        if( parsed_event_header && peek=='E' ) break;
    }

    // Check if all particles in last vertex were parsed
    /** @bug HepMC2 files produced with Pythia8 are known to have wrong
             information about number of particles in vertex. Hence '<' sign */
    if( is_parsing_successful && current_vertex_particles_parsed < current_vertex_particles_count ) {
        HEPMC3_ERROR( "ReaderAsciiHepMC2: not all particles parsed" )
        is_parsing_successful = false;
    }
    // Check if all vertices were parsed
    else if( is_parsing_successful && m_vertex_cache.size() != vertices_count ) {
        HEPMC3_ERROR( "ReaderAsciiHepMC2: not all vertices parsed" )
        is_parsing_successful = false;
    }

    if( !is_parsing_successful ) {
        HEPMC3_ERROR( "ReaderAsciiHepMC2: event parsing failed. Returning empty event" )
        HEPMC3_DEBUG( 1, "Parsing failed at line:" << std::endl << buf )
        evt.clear();
        m_isstream ? m_stream->clear(std::ios::badbit) : m_file.clear(std::ios::badbit);
        return 0;
    }

    // Restore production vertex pointers
    for(unsigned int i=0; i<m_particle_cache.size(); ++i) {
        if( !m_end_vertex_barcodes[i] ) continue;

        for(unsigned int j=0; j<m_vertex_cache.size(); ++j) {
            if( m_vertex_barcodes[j] == m_end_vertex_barcodes[i] ) {
                m_vertex_cache[j]->add_particle_in(m_particle_cache[i]);
                break;
            }
        }
    }

    // Remove vertices with no incoming particles or no outgoing particles
    for(unsigned int i=0; i<m_vertex_cache.size(); ++i) {
        if( m_vertex_cache[i]->particles_in().size() == 0 ) {
            HEPMC3_DEBUG( 30, "ReaderAsciiHepMC2::read_event - found a vertex without incoming particles: "<<m_vertex_cache[i]->id() );
//Sometimes the root vertex has no incoming particles.  Here we try to save the event.
            std::vector<GenParticlePtr> beams;
            for (auto p: m_vertex_cache[i]->particles_out()) if (p->status()==4 && !(p->end_vertex())) beams.push_back(p);
            for (auto p: beams)
            {
                m_vertex_cache[i]->add_particle_in(p);
                m_vertex_cache[i]->remove_particle_out(p);
                HEPMC3_DEBUG( 30, "ReaderAsciiHepMC2::read_event - moved particle with status=4 from the outgoing to the incoming particles of vertex: "<<m_vertex_cache[i]->id() );
            }
            if (beams.size()==0) m_vertex_cache[i] = nullptr;
        }
        else if( m_vertex_cache[i]->particles_out().size() == 0 ) {
            HEPMC3_DEBUG( 30, "ReaderAsciiHepMC2::read_event - found a vertex without outgouing particles: "<<m_vertex_cache[i]->id() );
            m_vertex_cache[i] = nullptr;
        }
    }

    // Reserve memory for the event
    evt.reserve( m_particle_cache.size(), m_vertex_cache.size() );

    // Add whole event tree in topological order
    evt.add_tree( m_particle_cache );

    for(unsigned int i=0; i<m_particle_cache.size(); ++i) {
        if(m_particle_cache_ghost[i]->attribute_names().size())
        {
            std::shared_ptr<DoubleAttribute> phi = m_particle_cache_ghost[i]->attribute<DoubleAttribute>("phi");
            if (phi) m_particle_cache[i]->add_attribute("phi",phi);
            std::shared_ptr<DoubleAttribute> theta = m_particle_cache_ghost[i]->attribute<DoubleAttribute>("theta");
            if (theta) m_particle_cache[i]->add_attribute("theta",theta);
            if (m_options.find("particle_flows_are_separated")!=m_options.end())
            {
                std::shared_ptr<IntAttribute> flow1 = m_particle_cache_ghost[i]->attribute<IntAttribute>("flow1");
                if (flow1) m_particle_cache[i]->add_attribute("flow1",flow1);
                std::shared_ptr<IntAttribute> flow2 = m_particle_cache_ghost[i]->attribute<IntAttribute>("flow2");
                if (flow2) m_particle_cache[i]->add_attribute("flow2",flow2);
                std::shared_ptr<IntAttribute> flow3 = m_particle_cache_ghost[i]->attribute<IntAttribute>("flow3");
                if (flow3) m_particle_cache[i]->add_attribute("flow3",flow3);
            }
            else
            {
                std::shared_ptr<VectorIntAttribute> flows = m_particle_cache_ghost[i]->attribute<VectorIntAttribute>("flows");
                if (flows)  m_particle_cache[i]->add_attribute("flows",flows);
            }
        }
    }

    for(unsigned int i=0; i<m_vertex_cache.size(); ++i)
        if(m_vertex_cache_ghost[i]->attribute_names().size())
        {
            for (size_t ii=0; ii<max_weights_size; ii++)
            {
                std::shared_ptr<DoubleAttribute> rs=m_vertex_cache_ghost[i]->attribute<DoubleAttribute>("weight"+std::to_string((long long unsigned int)ii));
                if (!rs) break;
                m_vertex_cache[i]->add_attribute("weight"+std::to_string((long long unsigned int)ii),rs);
            }
        }
    std::shared_ptr<IntAttribute> signal_process_vertex_barcode=evt.attribute<IntAttribute>("signal_process_vertex");
    if (signal_process_vertex_barcode) {
        int signal_process_vertex_barcode_value=signal_process_vertex_barcode->value();
        for(unsigned int i=0; i<m_vertex_cache.size(); ++i)
        {
            if (i>=m_vertex_barcodes.size()) continue;//this should not happen!
            if (signal_process_vertex_barcode_value!=m_vertex_barcodes.at(i)) continue;
            std::shared_ptr<IntAttribute> signal_process_vertex = std::make_shared<IntAttribute>(m_vertex_cache.at(i)->id());
            evt.add_attribute("signal_process_vertex",signal_process_vertex);
            break;
        }
    }
    m_particle_cache_ghost.clear();
    m_vertex_cache_ghost.clear();
    m_event_ghost->clear();
    return 1;
}

int ReaderAsciiHepMC2::parse_event_information(GenEvent &evt, const char *buf) {
    const char          *cursor             = buf;
    int                  event_no           = 0;
    int                  vertices_count     = 0;
    int                  random_states_size = 0;
    int                  weights_size       = 0;
    std::vector<long>    random_states(0);
    std::vector<double>  weights(0);

    // event number
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    event_no = atoi(cursor);
    evt.set_event_number(event_no);

    //mpi
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    std::shared_ptr<IntAttribute> mpi = std::make_shared<IntAttribute>(atoi(cursor));
    evt.add_attribute("mpi",mpi);

    //event scale
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    std::shared_ptr<DoubleAttribute> event_scale = std::make_shared<DoubleAttribute>(atof(cursor));
    evt.add_attribute("event_scale",event_scale);

    //alpha_qcd
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    std::shared_ptr<DoubleAttribute> alphaQCD = std::make_shared<DoubleAttribute>(atof(cursor));
    evt.add_attribute("alphaQCD",alphaQCD);

    //alpha_qed
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    std::shared_ptr<DoubleAttribute> alphaQED = std::make_shared<DoubleAttribute>(atof(cursor));
    evt.add_attribute("alphaQED",alphaQED);

    //signal_process_id
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    std::shared_ptr<IntAttribute> signal_process_id = std::make_shared<IntAttribute>(atoi(cursor));
    evt.add_attribute("signal_process_id",signal_process_id);

    //signal_process_vertex
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    std::shared_ptr<IntAttribute> signal_process_vertex = std::make_shared<IntAttribute>(atoi(cursor));
    evt.add_attribute("signal_process_vertex",signal_process_vertex);

    // num_vertices
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    vertices_count = atoi(cursor);

    // SKIPPED: beam 1
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;

    // SKIPPED: beam 2
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;

    //random states
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    random_states_size = atoi(cursor);
    random_states.resize(random_states_size);

    for ( int i = 0; i < random_states_size; ++i ) {
        if( !(cursor = strchr(cursor+1,' ')) ) return -1;
        random_states[i] = atoi(cursor);
    }
    if (m_options.find("event_random_states_are_separated")!=m_options.end())
    {
        evt.add_attribute("random_states",std::make_shared<VectorLongIntAttribute>(random_states));
    }
    else
    {
        for ( int i = 0; i < random_states_size; ++i )
            evt.add_attribute("random_states"+std::to_string((long long unsigned int)i),std::make_shared<IntAttribute>(random_states[i]));
    }
    // weights
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    weights_size = atoi(cursor);
    weights.resize(weights_size);

    for ( int i = 0; i < weights_size; ++i ) {
        if( !(cursor = strchr(cursor+1,' ')) ) return -1;
        weights[i] = atof(cursor);
    }

    evt.weights() = weights;

    HEPMC3_DEBUG( 10, "ReaderAsciiHepMC2: E: "<<event_no<<" ("<<vertices_count<<"V, "<<weights_size<<"W, "<<random_states_size<<"RS)" )

    return vertices_count;
}

bool ReaderAsciiHepMC2::parse_units(GenEvent &evt, const char *buf) {
    const char *cursor  = buf;

    // momentum
    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    ++cursor;
    Units::MomentumUnit momentum_unit = Units::momentum_unit(cursor);

    // length
    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    ++cursor;
    Units::LengthUnit length_unit = Units::length_unit(cursor);

    evt.set_units(momentum_unit,length_unit);

    HEPMC3_DEBUG( 10, "ReaderAsciiHepMC2: U: " << Units::name(evt.momentum_unit()) << " " << Units::name(evt.length_unit()) )

    return true;
}

int ReaderAsciiHepMC2::parse_vertex_information(const char *buf) {
    GenVertexPtr  data = std::make_shared<GenVertex>();
    GenVertexPtr  data_ghost = std::make_shared<GenVertex>();
    FourVector    position;
    const char   *cursor            = buf;
    int           barcode           = 0;
    int           num_particles_out = 0;
    int                  weights_size       = 0;
    std::vector<double>  weights(0);
    // barcode
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    barcode = atoi(cursor);

    // status
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    data->set_status( atoi(cursor) );

    // x
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    position.setX(atof(cursor));

    // y
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    position.setY(atof(cursor));

    // z
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    position.setZ(atof(cursor));

    // t
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    position.setT(atof(cursor));
    data->set_position( position );

    // SKIPPED: num_orphans_in
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;

    // num_particles_out
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    num_particles_out = atoi(cursor);

    //  weights

    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    weights_size = atoi(cursor);
    weights.resize(weights_size);

    for ( int i = 0; i < weights_size; ++i ) {
        if( !(cursor = strchr(cursor+1,' ')) ) return -1;
        weights[i] = atof(cursor);
    }



    // Add original vertex barcode to the cache
    m_vertex_cache.push_back( data );
    m_vertex_barcodes.push_back( barcode );

    m_event_ghost->add_vertex(data_ghost);
    if (m_options.find("vertex_weights_are_separated")!=m_options.end())
    {
        for ( int i = 0; i < weights_size; ++i )
            data_ghost->add_attribute("weight"+std::to_string((long long unsigned int)i),std::make_shared<DoubleAttribute>(weights[i]));
    }
    else
    {
        data_ghost->add_attribute("weights",std::make_shared<VectorDoubleAttribute>(weights));
        data_ghost->add_attribute("weights",std::make_shared<VectorDoubleAttribute>(weights));
    }
    m_vertex_cache_ghost.push_back( data_ghost );

    HEPMC3_DEBUG( 10, "ReaderAsciiHepMC2: V: "<<-(int)m_vertex_cache.size()<<" (old barcode"<<barcode<<") "<<num_particles_out<<" particles)" )

    return num_particles_out;
}

int ReaderAsciiHepMC2::parse_particle_information(const char *buf) {
    GenParticlePtr  data = std::make_shared<GenParticle>();
    GenParticlePtr  data_ghost = std::make_shared<GenParticle>();
    m_event_ghost->add_particle(data_ghost);
    FourVector      momentum;
    const char     *cursor  = buf;
    int             end_vtx = 0;

    /// @note barcode is ignored
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;

    // id
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    data->set_pid( atoi(cursor) );

    // px
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    momentum.setPx(atof(cursor));

    // py
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    momentum.setPy(atof(cursor));

    // pz
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    momentum.setPz(atof(cursor));

    // pe
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    momentum.setE(atof(cursor));
    data->set_momentum(momentum);

    // m
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    data->set_generated_mass( atof(cursor) );

    // status
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    data->set_status( atoi(cursor) );

    //theta
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    std::shared_ptr<DoubleAttribute> theta = std::make_shared<DoubleAttribute>(atof(cursor));
    if (theta->value()!=0.0) data_ghost->add_attribute("theta",theta);

    //phi
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    std::shared_ptr<DoubleAttribute> phi = std::make_shared<DoubleAttribute>(atof(cursor));
    if (phi->value()!=0.0) data_ghost->add_attribute("phi",phi);

    // end_vtx_code
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    end_vtx = atoi(cursor);

    //flow
    if( !(cursor = strchr(cursor+1,' ')) ) return -1;
    int flowsize=atoi(cursor);

    std::map<int,int> flows;
    for (int i=0; i<flowsize; i++)
    {
        if( !(cursor = strchr(cursor+1,' ')) ) return -1;
        int  flowindex=atoi(cursor);
        if( !(cursor = strchr(cursor+1,' ')) ) return -1;
        int flowvalue=atoi(cursor);
        flows[flowindex]=flowvalue;
    }
    if (m_options.find("particle_flows_are_separated")==m_options.end())
    {
        std::vector<int> vectorflows;
        for (auto f: flows) vectorflows.push_back(f.second);
        data_ghost->add_attribute("flows",std::make_shared<VectorIntAttribute>(vectorflows));
    }
    else
    {
        for (auto f: flows)   data_ghost->add_attribute("flow"+std::to_string((long long int)f.first),std::make_shared<IntAttribute>(f.second));
    }
// Set prod_vtx link
    if( end_vtx == m_vertex_barcodes.back() ) {
        m_vertex_cache.back()->add_particle_in(data);
        end_vtx = 0;
    }
    else {
        m_vertex_cache.back()->add_particle_out(data);
    }

    m_particle_cache.push_back( data );
    m_particle_cache_ghost.push_back( data_ghost );
    m_end_vertex_barcodes.push_back( end_vtx );

    HEPMC3_DEBUG( 10, "ReaderAsciiHepMC2: P: "<<m_particle_cache.size()<<" ( pid: "<<data->pid()<<") end vertex: "<<end_vtx )

    return 0;
}

bool ReaderAsciiHepMC2::parse_xs_info(GenEvent &evt, const char *buf) {
    const char *cursor  = buf;
    std::shared_ptr<GenCrossSection>  xs     = std::make_shared<GenCrossSection>();

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    double xs_val  = atof(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    double xs_err = atof(cursor);

    xs->set_cross_section( xs_val, xs_err);
    evt.add_attribute("GenCrossSection",xs);

    return true;
}

bool ReaderAsciiHepMC2::parse_weight_names(const char *buf) {
    const char     *cursor  = buf;
    const char     *cursor2 = buf;
    int             w_count = 0;
    std::vector<std::string>  w_names;

    // Ignore weight names if no GenRunInfo object
    if( !run_info() ) return true;

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    w_count = atoi(cursor);

    if( w_count <= 0 ) return false;

    w_names.resize(w_count);

    for( int i=0; i < w_count; ++i ) {
        // Find pair of '"' characters
        if( !(cursor  = strchr(cursor+1,'"')) ) return false;
        if( !(cursor2 = strchr(cursor+1,'"')) ) return false;

        // Strip cursor of leading '"' character
        ++cursor;

        w_names[i].assign(cursor, cursor2-cursor);

        cursor = cursor2;
    }

    run_info()->set_weight_names(w_names);

    return true;
}

bool ReaderAsciiHepMC2::parse_heavy_ion(GenEvent &evt, const char *buf) {
    std::shared_ptr<GenHeavyIon>  hi     = std::make_shared<GenHeavyIon>();
    const char              *cursor = buf;

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    hi->Ncoll_hard = atoi(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    hi->Npart_proj = atoi(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    hi->Npart_targ = atoi(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    hi->Ncoll = atoi(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    hi->spectator_neutrons = atoi(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    hi->spectator_protons = atoi(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    hi->N_Nwounded_collisions = atoi(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    hi->Nwounded_N_collisions = atoi(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    hi->Nwounded_Nwounded_collisions = atoi(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    hi->impact_parameter = atof(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    hi->event_plane_angle = atof(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    hi->eccentricity = atof(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    hi->sigma_inel_NN = atof(cursor);

    // Not in HepMC2:
    hi->centrality = 0.0;

    evt.add_attribute("GenHeavyIon",hi);

    return true;
}

bool ReaderAsciiHepMC2::parse_pdf_info(GenEvent &evt, const char *buf) {
    std::shared_ptr<GenPdfInfo>  pi     = std::make_shared<GenPdfInfo>();
    const char             *cursor = buf;

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    pi->parton_id[0] = atoi(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    pi->parton_id[1] = atoi(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    pi->x[0] = atof(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    pi->x[1] = atof(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    pi->scale = atof(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    pi->xf[0] = atof(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    pi->xf[1] = atof(cursor);

    //For compatibility with original HepMC2
    bool pdfids=true;
    if( !(cursor = strchr(cursor+1,' ')) ) pdfids=false;
    if(pdfids) pi->pdf_id[0] = atoi(cursor);
    else  pi->pdf_id[0] =0;

    if(pdfids) if( !(cursor = strchr(cursor+1,' ')) )  pdfids=false;
    if(pdfids) pi->pdf_id[1] = atoi(cursor);
    else  pi->pdf_id[1] =0;

    evt.add_attribute("GenPdfInfo",pi);

    return true;
}
bool ReaderAsciiHepMC2::failed() { return m_isstream ? (bool)m_stream->rdstate() :(bool)m_file.rdstate(); }

void ReaderAsciiHepMC2::close() {
    if (m_event_ghost) { m_event_ghost->clear(); delete m_event_ghost; m_event_ghost=nullptr;}
    if( !m_file.is_open() ) return;
    m_file.close();
}

} // namespace HepMC3
