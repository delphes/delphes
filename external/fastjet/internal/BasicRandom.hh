// Simple random number generator class taken from nlojet++.
// Some doxygen-style comments added by Gavin Salam, August 2006.
// $Id: BasicRandom.hh 1761 2010-09-16 10:43:18Z soyez $
//
//  Copyright (C) 2002 Zoltan Nagy
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#ifndef __FASTJET_BASICRANDOM_HH__
#define __FASTJET_BASICRANDOM_HH__ 1

//   Standard includes
#include <iostream>
#include <vector>
#include <cassert>
#include "fastjet/internal/base.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// \if internal_doc
/// @ingroup internal
/// \class BasicRandom
/// Base class for random number generator of a generic value type
/// \endif
template<typename _Tp> class BasicRandom {
public:
  typedef _Tp          value_type;
  typedef unsigned int size_type;
  typedef value_type*  pointer;
  
  //   give pseudo random numbers
  value_type operator() ();
  void       operator() (size_type, pointer);

  //   (re)initialize the random number generator
  void randomize(void *);

  //   minimum and maximum values
  static value_type min();
  static value_type max();

  //   print the informations about the generator to the stream
  void print_info(std::ostream& __os = std::cout);
};

//   default random generator
int __default_random_generator(int *__iseed); 


//   specializations

/// \if internal_doc
/// @ingroup internal
/// template specialization (int) for the BasicRandom template class. 
/// \endif
template<>
class BasicRandom<int>
{
public:
  typedef int          value_type;
  typedef unsigned int size_type;
  typedef value_type*  pointer;
  
  // constructors
  explicit BasicRandom(int __s1 = 12345, int __s2 = 67890) {
    _M_iseed[0] = __s1;
    _M_iseed[1] = __s2;
  }
    
  //   give pseudo random numbers
  value_type operator() () { 
    return __default_random_generator(_M_iseed);
  }
  
  void operator() (size_type __n, pointer __res) {
    for(size_type __i = 0; __i < __n; __i++) 
      __res[__i] = __default_random_generator(_M_iseed);
  }

  //   (re)initialize the random number generator
  void randomize(void *__iseed) {
    int *__new_seed = (int*) __iseed;
    _M_iseed[0] = __new_seed[0];
    _M_iseed[1] = __new_seed[1];
  }

  void set_status(const std::vector<int> & __iseed) {
    assert(__iseed.size() >= 2);
    _M_iseed[0] = __iseed[0];
    _M_iseed[1] = __iseed[1];
  }

  void get_status(std::vector<int> & __iseed) {
    __iseed.resize(2);
    __iseed[0] = _M_iseed[0];
    __iseed[1] = _M_iseed[1];
  }
  
  //   minimum and maximum values
  inline static value_type min() { return 0;}
  inline static value_type max() { return 2147483647;}

  //   print the informations about the generator to the stream
  void print_info(std::ostream& __os = std::cout) {
    __os<<"BasicRandom<int> : "<<_M_iseed[0]<<", "<<_M_iseed[1]<<std::endl;
  }
  
private:
  int _M_iseed[2];
};
  

/// \if internal_doc
/// @ingroup internal
/// template specialization (double) for the BasicRandom template class. 
/// \endif
template<> class BasicRandom<double> {
public:
  typedef double       value_type;
  typedef unsigned int size_type;
  typedef value_type*  pointer;
  
  /// constructor that takes two integers to specify the seed
  explicit BasicRandom(int __s1 = 12345, int __s2 = 67890) {
    _M_iseed[0] = __s1;
    _M_iseed[1] = __s2;
  }
    
  /// return a single pseudorandom double number, in the range 0.0 to 1.0
  /// (not sure whether this range is open or closed)
  value_type operator() () { 
    return 4.6566128752457969241e-10*__default_random_generator(_M_iseed);
  }
  
  /// given a pointer __res to the beginning of an array, fill that array
  /// with __n random numbers
  void operator() (size_type __n, pointer __res) {
    for(size_type __i = 0; __i < __n; __i++) 
      __res[__i] = this -> operator()(); 
  }

  ///  (re)initialize the random number generator from an array of seeds
  void randomize(void *__iseed) {
    int *__new_seed = (int*) __iseed;
    _M_iseed[0] = __new_seed[0];
    _M_iseed[1] = __new_seed[1];
  }
  
  void set_status(const std::vector<int> & __iseed) {
    assert(__iseed.size() >= 2);
    _M_iseed[0] = __iseed[0];
    _M_iseed[1] = __iseed[1];
  }

  void get_status(std::vector<int> & __iseed) {
    __iseed.resize(2);
    __iseed[0] = _M_iseed[0];
    __iseed[1] = _M_iseed[1];
  }
  
  /// minimum value returned by the generator
  inline static value_type min() { return 0.0;}
  /// maximum value returned by the generator
  inline static value_type max() { return 1.0;}

  ///  print information about the generator to the stream
  void print_info(std::ostream& __os = std::cout) {
    __os<<"BasicRandom<double> : "<<_M_iseed[0]<<", "<<_M_iseed[1]<<std::endl;
  }
  
private:
  int _M_iseed[2];
};
  
//   globally defined random number generator
extern BasicRandom<int>     _G_random_int;
extern BasicRandom<double>  _G_random_double;


FASTJET_END_NAMESPACE

#endif // __FASTJET_BASICRANDOM_HH__

