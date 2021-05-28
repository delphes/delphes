// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
///
/// @file ReaderAscii.cc
/// @brief Implementation of \b class ReaderAscii
///
#include "HepMC3/ReaderAscii.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Units.h"
#include <cstring>
#include <sstream>

namespace HepMC3 {


ReaderAscii::ReaderAscii(const std::string &filename)
    : m_file(filename), m_stream(0), m_isstream(false)
{
    if( !m_file.is_open() ) {
        HEPMC3_ERROR( "ReaderAscii: could not open input file: "<<filename )
    }
    set_run_info(std::make_shared<GenRunInfo>());
}


// Ctor for reading from stdin
ReaderAscii::ReaderAscii(std::istream & stream)
    : m_stream(&stream), m_isstream(true)
{
    if( !m_stream->good() ) {
        HEPMC3_ERROR( "ReaderAscii: could not open input stream " )
    }
    set_run_info(std::make_shared<GenRunInfo>());
}



ReaderAscii::~ReaderAscii() { if (!m_isstream) close(); }

bool ReaderAscii::skip(const int n)
{
    const size_t       max_buffer_size=512*512;
    char               buf[max_buffer_size];
    int nn=n;
    while(!failed()) {
        char  peek;
        if ( (!m_file.is_open()) && (!m_isstream) ) return false;
        m_isstream ? peek = m_stream->peek() : peek = m_file.peek();
        if( peek=='E' ) nn--;
        if (nn<0) return true;
        m_isstream ? m_stream->getline(buf,max_buffer_size) : m_file.getline(buf,max_buffer_size);
    }
    return true;
}


bool ReaderAscii::read_event(GenEvent &evt) {
    if ( (!m_file.is_open()) && (!m_isstream) ) return false;

    char               peek;
    const size_t       max_buffer_size=512*512;
    char               buf[max_buffer_size];
    bool               parsed_event_header    = false;
    bool               is_parsing_successful  = true;
    std::pair<int,int> vertices_and_particles(0,0);

    evt.clear();
    evt.set_run_info(run_info());
    m_forward_daughters.clear();
    m_forward_mothers.clear();
    //
    // Parse event, vertex and particle information
    //
    while(!failed()) {

        m_isstream ? m_stream->getline(buf,max_buffer_size) : m_file.getline(buf,max_buffer_size);

        if( strlen(buf) == 0 ) continue;

        // Check for ReaderAscii header/footer
        if( strncmp(buf,"HepMC",5) == 0 ) {
            if( strncmp(buf,"HepMC::Version",14) != 0 && strncmp(buf,"HepMC::Asciiv3",14)!=0 )
            {
                HEPMC3_WARNING( "ReaderAscii: found unsupported expression in header. Will close the input." )
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
            vertices_and_particles = parse_event_information(evt,buf);
            if (vertices_and_particles.second < 0) {
                is_parsing_successful = false;
            } else {
                is_parsing_successful = true;
                parsed_event_header   = true;
            }
            break;
        case 'V':
            is_parsing_successful = parse_vertex_information(evt,buf);
            break;
        case 'P':
            is_parsing_successful = parse_particle_information(evt,buf);
            break;
        case 'W':
            if ( parsed_event_header )
                is_parsing_successful = parse_weight_values(evt,buf);
            else
                is_parsing_successful = parse_weight_names(buf);
            break;
        case 'U':
            is_parsing_successful = parse_units(evt,buf);
            break;
        case 'T':
            is_parsing_successful = parse_tool(buf);
            break;
        case 'A':
            if ( parsed_event_header )
                is_parsing_successful = parse_attribute(evt,buf);
            else
                is_parsing_successful = parse_run_attribute(buf);
            break;
        default:
            HEPMC3_WARNING( "ReaderAscii: skipping unrecognised prefix: " << buf[0] )
            is_parsing_successful = true;
            break;
        }

        if( !is_parsing_successful ) break;

        // Check for next event
        m_isstream ? peek = m_stream->peek() : peek = m_file.peek();
        if( parsed_event_header && peek=='E' ) break;
    }


    // Check if all particles and vertices were parsed
    if ((int)evt.particles().size() > vertices_and_particles.second ) {
        HEPMC3_ERROR( "ReaderAscii: too many particles were parsed" )
        printf("%zu  vs  %i expected\n",evt.particles().size(),vertices_and_particles.second );
        is_parsing_successful = false;
    }
    if ((int)evt.particles().size() < vertices_and_particles.second ) {
        HEPMC3_ERROR( "ReaderAscii: too few  particles were parsed" )
        printf("%zu  vs  %i expected\n",evt.particles().size(),vertices_and_particles.second );
        is_parsing_successful = false;
    }

    if ((int)evt.vertices().size()  > vertices_and_particles.first) {
        HEPMC3_ERROR( "ReaderAscii: too many vertices were parsed" )
        printf("%zu  vs  %i expected\n",evt.vertices().size(),vertices_and_particles.first );
        is_parsing_successful =  false;
    }

    if ((int)evt.vertices().size()  < vertices_and_particles.first) {
        HEPMC3_ERROR( "ReaderAscii: too few vertices were parsed" )
        printf("%zu  vs  %i expected\n",evt.vertices().size(),vertices_and_particles.first );
        is_parsing_successful =  false;
    }
    // Check if there were HEPMC3_ERRORs during parsing
    if( !is_parsing_successful ) {
        HEPMC3_ERROR( "ReaderAscii: event parsing failed. Returning empty event" )
        HEPMC3_DEBUG( 1, "Parsing failed at line:" << std::endl << buf )

        evt.clear();
        m_isstream ? m_stream->clear(std::ios::badbit) : m_file.clear(std::ios::badbit);

        return false;
    }
    for ( auto p : m_forward_daughters )
        for (auto v: evt.vertices())
            if (p.second==v->id())
                v->add_particle_out(p.first);
    for ( auto v : m_forward_mothers )  for ( auto idpm : v.second )  v.first->add_particle_in(evt.particles()[idpm-1]);

    /* restore ids of vertices using a bank of available ids*/
    std::vector<int> all_ids;
    std::vector<int> filled_ids;
    std::vector<int> diff;
    for (auto v: evt.vertices()) if (v->id()!=0) filled_ids.push_back(v->id());
    for (int i=-((long)evt.vertices().size()); i<0; i++) all_ids.push_back(i);
    std::sort(all_ids.begin(),all_ids.end());
    std::sort(filled_ids.begin(),filled_ids.end());
    //The bank of available ids is created as a difference between all range of ids and the set of used ids
    std::set_difference(all_ids.begin(), all_ids.end(), filled_ids.begin(), filled_ids.end(), std::inserter(diff, diff.begin()));
    auto it= diff.rbegin();
    //Set available ids to vertices sequentially.
    for (auto v: evt.vertices()) if (v->id()==0) { v->set_id(*it); it++;}

    return true;
}


std::pair<int,int> ReaderAscii::parse_event_information(GenEvent &evt, const char *buf) {
    static const std::pair<int,int>  err(-1,-1);
    std::pair<int,int>               ret(-1,-1);
    const char                 *cursor   = buf;
    int                         event_no = 0;
    FourVector                  position;

    // event number
    if( !(cursor = strchr(cursor+1,' ')) ) return err;
    event_no = atoi(cursor);
    evt.set_event_number(event_no);

    // num_vertices
    if( !(cursor = strchr(cursor+1,' ')) ) return err;
    ret.first = atoi(cursor);

    // num_particles
    if( !(cursor = strchr(cursor+1,' ')) ) return err;
    ret.second = atoi(cursor);

    // check if there is position information
    if( (cursor = strchr(cursor+1,'@')) ) {

        // x
        if( !(cursor = strchr(cursor+1,' ')) ) return err;
        position.setX(atof(cursor));

        // y
        if( !(cursor = strchr(cursor+1,' ')) ) return err;
        position.setY(atof(cursor));

        // z
        if( !(cursor = strchr(cursor+1,' ')) ) return err;
        position.setZ(atof(cursor));

        // t
        if( !(cursor = strchr(cursor+1,' ')) ) return err;
        position.setT(atof(cursor));
        evt.shift_position_to( position );
    }

    HEPMC3_DEBUG( 10, "ReaderAscii: E: "<<event_no<<" ("<<ret.first<<"V, "<<ret.second<<"P)" )

    return ret;
}


bool ReaderAscii::parse_weight_values(GenEvent &evt, const char *buf) {

    std::istringstream iss(buf + 1);
    std::vector<double> wts;
    double w;
    while ( iss >> w ) wts.push_back(w);
    if ( run_info() && run_info()->weight_names().size()
            && run_info()->weight_names().size() != wts.size() )
        throw std::logic_error("ReaderAscii::parse_weight_values: "
                               "The number of weights ("+std::to_string((long long int)(wts.size()))+") does not match "
                               "the  number weight names("+std::to_string((long long int)(run_info()->weight_names().size()))+") in the GenRunInfo object");
    evt.weights() = wts;

    return true;
}


bool ReaderAscii::parse_units(GenEvent &evt, const char *buf) {
    const char *cursor = buf;

    // momentum
    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    ++cursor;
    Units::MomentumUnit momentum_unit = Units::momentum_unit(cursor);

    // length
    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    ++cursor;
    Units::LengthUnit length_unit = Units::length_unit(cursor);

    evt.set_units(momentum_unit,length_unit);

    HEPMC3_DEBUG( 10, "ReaderAscii: U: " << Units::name(evt.momentum_unit()) << " " << Units::name(evt.length_unit()) )

    return true;
}


bool ReaderAscii::parse_vertex_information(GenEvent &evt, const char *buf) {
    GenVertexPtr  data = std::make_shared<GenVertex>();
    FourVector    position;
    const char   *cursor          = buf;
    const char   *cursor2         = nullptr;
    int           id              = 0;
    int           highest_id      = evt.particles().size();

    // id
    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    id = atoi(cursor);

    // status
    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    data->set_status( atoi(cursor) );

    // skip to the list of particles
    if( !(cursor = strchr(cursor+1,'[')) ) return false;

    while(true) {
        ++cursor;             // skip the '[' or ',' character
        cursor2     = cursor; // save cursor position
        int  particle_in = atoi(cursor);

        // add incoming particle to the vertex
        if( particle_in > 0) {
//Particles are always ordered, so id==position in event.
            if (particle_in <= highest_id)
                data->add_particle_in( evt.particles()[particle_in-1] );
//If the particle has not been red yet, we store its id to add the particle later.
            else m_forward_mothers[data].insert(particle_in);
        }

        // check for next particle or end of particle list
        if( !(cursor = strchr(cursor+1,',')) ) {
            if( !(cursor = strchr(cursor2+1,']')) ) return false;
            break;
        }
    }

    // check if there is position information
    if( (cursor = strchr(cursor+1,'@')) ) {

        // x
        if( !(cursor = strchr(cursor+1,' ')) ) return false;
        position.setX(atof(cursor));

        // y
        if( !(cursor = strchr(cursor+1,' ')) ) return false;
        position.setY(atof(cursor));

        // z
        if( !(cursor = strchr(cursor+1,' ')) ) return false;
        position.setZ(atof(cursor));

        // t
        if( !(cursor = strchr(cursor+1,' ')) ) return false;
        position.setT(atof(cursor));
        data->set_position( position );

    }

    HEPMC3_DEBUG( 10, "ReaderAscii: V: "<<id<<" with "<<data->particles_in().size()<<" particles)" )

    evt.add_vertex(data);
//Restore vertex id, as it is used to build connections inside event.
    data->set_id(id);

    return true;
}


bool ReaderAscii::parse_particle_information(GenEvent &evt, const char *buf) {
    GenParticlePtr  data = std::make_shared<GenParticle>();
    FourVector      momentum;
    const char     *cursor  = buf;
    int             mother_id = 0;

    // verify id
    if( !(cursor = strchr(cursor+1,' ')) ) return false;

    if( atoi(cursor) != (int)evt.particles().size() + 1 ) {
        /// @todo Should be an exception
        HEPMC3_ERROR( "ReaderAscii: particle ID mismatch" )
        return false;
    }

    // mother id
    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    mother_id = atoi(cursor);

    // Parent object is a particle. Particleas are always ordered id==position in event.
    if( mother_id > 0 && mother_id <= (int)evt.particles().size() ) {

        GenParticlePtr mother = evt.particles()[ mother_id-1 ];
        GenVertexPtr   vertex = mother->end_vertex();

        // create new vertex if needed
        if( !vertex ) {
            vertex = std::make_shared<GenVertex>();
            vertex->add_particle_in(mother);
        }

        vertex->add_particle_out(data);
        evt.add_vertex(vertex);
//ID of this vertex is not explicitely set in the input. We set it to zero to prevent overlap with other ids. It will be restored later.
        vertex->set_id(0);
    }
    // Parent object is vertex
    else if( mother_id < 0 )
    {
        //Vertices are not always ordered, e.g. when one reads HepMC2 event, so we check their ids.
        bool found=false;
        for (auto v: evt.vertices()) if (v->id()==mother_id) {v->add_particle_out(data); found=true; break; }
        if (!found)
        {
//This should happen  in case of unordered event.
//      WARNING("ReaderAscii: Unordered event, id of mother vertex  is out of range of known ids:   " <<mother_id<<" evt.vertices().size()="<<evt.vertices().size() )
//Save the mother id to reconnect later.
            m_forward_daughters[data]=mother_id;
        }
    }

    // pdg id
    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    data->set_pid( atoi(cursor) );

    // px
    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    momentum.setPx(atof(cursor));

    // py
    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    momentum.setPy(atof(cursor));

    // pz
    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    momentum.setPz(atof(cursor));

    // pe
    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    momentum.setE(atof(cursor));
    data->set_momentum(momentum);

    // m
    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    data->set_generated_mass( atof(cursor) );

    // status
    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    data->set_status( atoi(cursor) );

    evt.add_particle(data);

    HEPMC3_DEBUG( 10, "ReaderAscii: P: "<<data->id()<<" ( mother: "<<mother_id<<", pid: "<<data->pid()<<")" )

    return true;
}


bool ReaderAscii::parse_attribute(GenEvent &evt, const char *buf) {
    const char     *cursor  = buf;
    const char     *cursor2 = buf;
    char            name[512];
    int             id = 0;

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    id = atoi(cursor);

    if( !(cursor  = strchr(cursor+1,' ')) ) return false;
    ++cursor;

    if( !(cursor2 = strchr(cursor,' ')) ) return false;
    snprintf(name, 512,"%.*s",(int)(cursor2-cursor), cursor);

    cursor = cursor2+1;

    std::shared_ptr<Attribute> att =
        std::make_shared<StringAttribute>( StringAttribute(unescape(cursor)) );

    evt.add_attribute(std::string(name), att, id);

    return true;
}

bool ReaderAscii::parse_run_attribute(const char *buf) {
    const char     *cursor  = buf;
    const char     *cursor2 = buf;
    char            name[512];

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    ++cursor;

    if( !(cursor2 = strchr(cursor,' ')) ) return false;
    snprintf(name, 512,"%.*s", (int)(cursor2-cursor), cursor);

    cursor = cursor2+1;

    std::shared_ptr<StringAttribute> att =
        std::make_shared<StringAttribute>( StringAttribute(unescape(cursor)) );

    run_info()->add_attribute(std::string(name), att);

    return true;

}


bool ReaderAscii::parse_weight_names(const char *buf) {
    const char     *cursor  = buf;

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    ++cursor;

    std::istringstream iss(unescape(cursor));
    std::vector<std::string> names;
    std::string name;
    while ( iss >> name ) names.push_back(name);

    run_info()->set_weight_names(names);

    return true;

}

bool ReaderAscii::parse_tool(const char *buf) {
    const char     *cursor  = buf;

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    ++cursor;
    std::string line = unescape(cursor);
    GenRunInfo::ToolInfo tool;
    std::string::size_type pos = line.find("\n");
    tool.name = line.substr(0, pos);
    line = line.substr(pos + 1);
    pos = line.find("\n");
    tool.version = line.substr(0, pos);
    tool.description = line.substr(pos + 1);
    run_info()->tools().push_back(tool);

    return true;

}


std::string ReaderAscii::unescape(const std::string& s) {
    std::string ret;
    ret.reserve(s.length());
    for ( std::string::const_iterator it = s.begin(); it != s.end(); ++it ) {
        if ( *it == '\\' ) {
            ++it;
            if ( *it == '|' )
                ret += '\n';
            else
                ret += *it;
        } else
            ret += *it;
    }

    return ret;
}

bool ReaderAscii::failed() { return m_isstream ? (bool)m_stream->rdstate() :(bool)m_file.rdstate(); }

void ReaderAscii::close() {
    if( !m_file.is_open()) return;
    m_file.close();
}


} // namespace HepMC3
