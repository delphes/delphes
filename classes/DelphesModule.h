/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DelphesModule_h
#define DelphesModule_h

/** \class DelphesModule
 *
 *  Base class for all Delphes modules
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModuleFactory.h"
#include "classes/DelphesParameters.h"

class Candidate;

class DelphesFactory;

class DelphesModule
{
public:
  explicit DelphesModule(const DelphesParameters &);
  virtual ~DelphesModule() = default;

  virtual void Init() {}
  virtual void Process() {}
  virtual void Finish() {}

  std::shared_ptr<std::vector<Candidate *> > ImportArray(std::string_view name);
  std::shared_ptr<std::vector<Candidate *> > ExportArray(std::string_view name);

  void SetName(std::string_view moduleName) { fName = moduleName; }
  const std::string &GetName() const { return fName; }

  virtual void SetFactory(DelphesFactory *factory) { fFactory = factory; }
  virtual DelphesFactory *GetFactory() const;

  virtual bool IsReader() const { return false; }
  virtual bool IsWriter() const { return false; }

protected:
  template <typename T>
  T Steer(std::string_view keyName, T defaultValue = T{}) const
  {
    return fModuleParams.Get<T>(std::string{keyName}, defaultValue);
  }
  std::string Steer(std::string_view keyName, std::string_view defaultValue = {}) const
  {
    return fModuleParams.Get<std::string>(std::string{keyName}, std::string{defaultValue});
  }

private:
  std::string fName;
  DelphesFactory *fFactory{nullptr};

  const DelphesParameters fModuleParams;
};

/// Add a processing module to the list of handled modules
#define REGISTER_MODULE(name, obj)                                                           \
  struct BUILDER_NAME(obj)                                                                   \
  {                                                                                          \
    BUILDER_NAME(obj)() { DelphesProcessingModuleFactory::Get().RegisterModule<obj>(name); } \
  };                                                                                         \
  static const BUILDER_NAME(obj) gDelphesModule##obj;                                        \
  static_assert(true, "")

/// A documentation generator factory
DEFINE_FACTORY(DelphesProcessingModuleFactory, DelphesModule, "Processing modules factory");

#endif /* DelphesModule_h */
