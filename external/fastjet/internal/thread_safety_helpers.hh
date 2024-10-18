//FJSTARTHEADER
// $Id$
//
// Copyright (c) 2014-2024, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development. They are described in the original FastJet paper,
//  hep-ph/0512210 and in the manual, arXiv:1111.6097. If you use
//  FastJet as part of work towards a scientific publication, please
//  quote the version you use and include a citation to the manual and
//  optionally also to hep-ph/0512210.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//FJENDHEADER

#ifndef __FASTJET_THREAD_SAFETY_HELPERS_HH__
#define __FASTJET_THREAD_SAFETY_HELPERS_HH__

/// The code in this file is supposed to help writing code that will
/// automatically provide thread-safe features when available and come
/// revert back to "old/standard" C++ if thread-safety is not switched
/// on
///
///\TODO fix doxygen comments (declare things as internal; make sure
/// doxygen doc is not duplicate --- if necessary, keep only doxygen
/// comments in the thread-safe versions)

#include "fastjet/internal/base.hh"
#include "fastjet/config.h"
#include <limits>

#ifdef FASTJET_HAVE_LIMITED_THREAD_SAFETY

// introduces a few tools in CXX11 that we'll use in some FJ classes
#include <atomic>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace thread_safety_helpers{

  //----------------------------------------------------------------------
  /// \if internal_doc
  /// \class AtomicCounter
  ///
  /// provides a thread-safe counter which can only step one unit at a time
  /// and has overflow protection
  /// \endif
  template<typename T>
  class AtomicCounter{
  public:
    /// default ctor
    AtomicCounter() : _count{0}{}

    /// ctor with initialisation
    AtomicCounter(const T &count) : _count{count}{}

    /// copy ctor (works around the deleted copy in atomic, see
    /// e.g. http://stackoverflow.com/questions/19883092/error-implicitly-deleted-because-the-default-definition-would-be-ill-formed-ve)
    AtomicCounter(const AtomicCounter &other) : _count{other._count.load()}{}

    /// for a more friendly usage, overload the type cast to the
    /// base template type
    operator T() const{ return _count.load();}

    /// get the count
    T get() const{ return _count.load();}
    
    /// set the counter to a given value
    void set(const T new_value){
      _count.store(new_value);
    }

    /// step the counter and return the count just before it was stepped
    ///
    /// Q: can we declare this as T && ...?
    T step(){
      // another thread could be upadting this at the same time, so extra
      // care is needed.
      //
      // Recall that the compare_exchange_strong will return true if the
      // exchange has been done. Otherwise, it means that the count
      // changed in the meantime, so we try again. Also, since when it
      // "fails" compare_exchange_strong loads the count of *this in
      // expected, count does not need to be re-read in the loop!
      //
      // Note that at the end of this procedure, count will countain the
      // number of times this warning occured just before this
      // occurence. It can thus be used to see if it needs to be printed
      // out
      //
      // Note also that compared to the apparently simpler fetch_add,
      // this method also avoids overflows
      T count = _count;
      while (_count < std::numeric_limits<T>::max()
             && !(_count.compare_exchange_strong(count, count+1)));
      return count;
    }

    /// override the ++ operator
    /// prefix version
    inline T operator++(){
      return step()+1;
    }

    /// override the ++ operator
    /// postfix version
    inline T operator++(int){
      return step();
    }

  private:
    std::atomic<T> _count;  ///< the actual count
  };
  
  //----------------------------------------------------------------------
  /// \if internal_doc
  /// \class FirstTimeTrue
  /// provides an object wich will return "true" the first time () is
  /// called and false afterwards
  /// \endif
  class FirstTimeTrue{
  public:
    FirstTimeTrue(): _first_time{true}{}
    // explicit copy ctor (this class contains atimoc vars)
    FirstTimeTrue(const FirstTimeTrue &other) : _first_time{other._first_time.load()}{}
    bool operator()(){
      // Thread-safety note:
      //   the construct
      //      if (!_first_time) {return;}
      //      _first_time = false;
      //   is dangerous because the test can be passed by a second thread
      //   before the first one has set it to false. Use atomic exchange
      //   to handle this better
      bool expected = true;
      // this behaves as follows: if we have the expected value (true),
      // set _first_time to the desired (false) and return
      // true. Otherwise, do nothing and return false
      //
      // Note that since we are not using the "expected" value
      // afterwards, we can use a relaxed memory ordering if the next
      // call returns false
      return _first_time.compare_exchange_strong(expected, false,
                                                 std::memory_order_seq_cst,
                                                 std::memory_order_relaxed);
    }
  private:
    std::atomic<bool> _first_time;
  };

} // namespace thread_safety_helpers

FASTJET_END_NAMESPACE

#else  // FJ wo thread-safety features

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace thread_safety_helpers{
  //----------------------------------------------------------------------
  /// \class AtomicCounter
  ///
  /// (would) provides a thread-safe counter (with CXX11 features)
  template<typename T>
  class AtomicCounter{
  public:
    /// default ctor
    AtomicCounter() : _count(0){}

    /// ctor with initialisation
    AtomicCounter(const T &count) : _count(count){}

    /// copy ctor
    AtomicCounter(const AtomicCounter &other) : _count(other._count){}

    /// for a more friendly usage, overload the type cast
    ///
    /// This will (likely) allow a transparent usage w or wo C++11
    /// features enabled
    operator T() const{ return _count;}
       
    /// get the count
    T get() const{ return _count;}

    /// set the counter to a given value
    void set(const T new_value){
      _count = new_value;
    }

    /// step the counter and return the value just before it was stepped
    T step(){
      unsigned int count = _count;
      if (_count < std::numeric_limits<T>::max()){ _count++; }
      return count;
    }

    /// override the ++ operator
    /// prefix version
    inline T operator++(){
      return step()+1;
    }

    /// override the ++ operator
    /// postfix version
    inline T operator++(int){
      return step();
    }
    
  private:
    T _count;  ///< the actual value
  };

  //----------------------------------------------------------------------
  /// \class FirstTimeTrue
  /// provides an object wich will return "true" the first time () is
  /// called and false afterwards
  class FirstTimeTrue{
  public:
    FirstTimeTrue(): _first_time(true){}
    bool operator()(){
      if (!_first_time) {return false;}
      _first_time = false;
      return true;
    }
  private:
    bool _first_time;
  };
} // namespace thread_safety_helpers

FASTJET_END_NAMESPACE



#endif // FASTJET_HAVE_LIMITED_THREAD_SAFETY

#endif // __FASTJET_THREAD_SAFETY_HELPERS_HH__
