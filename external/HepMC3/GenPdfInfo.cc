// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file GenPdfInfo.cc
 *  @brief Implementation of \b class GenPdfInfo
 *
 */
#include "HepMC3/GenPdfInfo.h"
#include <cstring> // memcmp
#include <cstdlib> // atoi
#include <cstdio> // sprintf

namespace HepMC3 {

bool GenPdfInfo::from_string(const std::string &att) {
    const char *cursor = att.data();

    parton_id[0] = atoi(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    parton_id[1] = atoi(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    x[0] = atof(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    x[1] = atof(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    scale = atof(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    xf[0] = atof(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    xf[1] = atof(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    pdf_id[0] = atoi(cursor);

    if( !(cursor = strchr(cursor+1,' ')) ) return false;
    pdf_id[1] = atoi(cursor);

    return true;
}

bool GenPdfInfo::to_string(std::string &att) const {
    char buf[255];//Note: the format is fixed, so no reason for complicatied tratment

    snprintf(buf,255,"%i %i %.8e %.8e %.8e %.8e %.8e %i %i",
             parton_id[0],
             parton_id[1],
             x[0],
             x[1],
             scale,
             xf[0],
             xf[1],
             pdf_id[0],
             pdf_id[1]);

    att = buf;

    return true;
}

void GenPdfInfo::set(const int& parton_id1, const int& parton_id2, const double& x1, const double& x2,
                     const double& scale_in, const double& xf1,const double& xf2,
                     const int& pdf_id1, const int& pdf_id2) {
    parton_id[0] = parton_id1;
    parton_id[1] = parton_id2;
    x[0]         = x1;
    x[1]         = x2;
    scale        = scale_in;
    xf[0]        = xf1;
    xf[1]        = xf2;
    pdf_id[0]    = pdf_id1;
    pdf_id[1]    = pdf_id2;
}

bool GenPdfInfo::operator==( const GenPdfInfo& a ) const {
    return ( memcmp( (void*)this, (void*)&a, sizeof(class GenPdfInfo) ) == 0 );
}

bool GenPdfInfo::operator!=( const GenPdfInfo& a ) const {
    return !( a == *this );
}

bool GenPdfInfo::is_valid() const
{
    if( parton_id[0] != 0 ) return true;
    if( parton_id[1] != 0 ) return true;
    if( x[0]         != 0 ) return true;
    if( x[1]         != 0 ) return true;
    if( scale        != 0 ) return true;
    if( xf[0]        != 0 ) return true;
    if( xf[1]        != 0 ) return true;
    if( pdf_id[0]    != 0 ) return true;
    if( pdf_id[1]    != 0 ) return true;

    return false;
}

} // namespace HepMC3
