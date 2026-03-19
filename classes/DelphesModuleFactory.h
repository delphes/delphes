/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2026  Universite catholique de Louvain (UCL), Belgium
 *                           AGH University of Krakow, Poland
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

#ifndef classes_DelphesModuleFactory_h
#define classes_DelphesModuleFactory_h

#include <algorithm>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

/// Name of the object Builder
#define BUILDER_NAME(obj) obj##Builder
/// Define a new factory instance for the definition of modules
#define DEFINE_FACTORY(name, obj_type, type)        \
  class name: public DelphesModuleFactory<obj_type> \
  {                                                 \
  public:                                           \
    inline static name &Get()                       \
    {                                               \
      static name instance;                         \
      return instance;                              \
    }                                               \
                                                    \
  private:                                          \
    explicit name() : DelphesModuleFactory(type) {} \
  };                                                \
  static_assert(true, "")

template <typename T>
class DelphesModuleFactory
{
public:
  DelphesModuleFactory(const DelphesModuleFactory &) = delete; ///< Disabled copy constructor
  virtual ~DelphesModuleFactory() = default; ///< Default destructor
  void operator=(const DelphesModuleFactory &) = delete; ///< Disabled assignment operator

  static DelphesModuleFactory &Get()
  {
    static DelphesModuleFactory instance;
    return instance;
  }
  /// Register a named module in the database
  template <typename U>
  void RegisterModule(const std::string &name)
  {
    static_assert(std::is_base_of_v<T, U>,
      "\n\n  *** Failed to register a module with improper inheritance into the factory. ***\n");
    if(Has(name))
      ThrowDuplicateError(name);
    fBuildersMap.insert(std::make_pair(name, &BuildModule<U>));
  }
  /// Build one instance of a named module
  std::unique_ptr<T> Build(const std::string &name) const
  {
    if(fBuildersMap.count(name) > 0)
      return fBuildersMap.at(name)();
    ThrowBuildError(name);
    throw;
  }

  typedef std::unique_ptr<T> (*Builder)(); ///< constructor type for a module definition
  bool Has(const std::string &name) const { return fBuildersMap.count(name) > 0; } ///< is a named module registered?
  const std::string &Type() const { return fType; } ///< module types in factory
  /// List of module names registered in the factory
  std::vector<std::string> Modules() const
  {
    std::vector<std::string> modules;
    for(const auto &mod : fBuildersMap)
      modules.emplace_back(mod.first);
    std::sort(modules.begin(), modules.end());
    return modules;
  }

protected:
  explicit DelphesModuleFactory(std::string type) : fType(std::move(type)) {}

private:
  void ThrowDuplicateError(std::string_view name) const
  {
    std::ostringstream message;
    message << "\n\n  *** Factory of " << fType << " modules detected a duplicate " << fType << " registration for name '" << name << "'! ***\n";
    throw std::invalid_argument(message.str());
  }
  void ThrowBuildError(std::string_view name) const
  {
    std::ostringstream message;
    message << "\n\n  *** Failed to build a " << fType << " module with name '" << std::string{name} + "' ***\n\n"
            << "     List of " + fType + " modules handled:\n";
    for(const auto &module_name : Modules())
      message << "     * " << module_name << "\n";
    throw std::invalid_argument(message.str());
  }

  /// Build a module with its parameters set
  template <typename U>
  static std::unique_ptr<T> BuildModule()
  {
    return std::make_unique<U>();
  }
  std::unordered_map<std::string, Builder> fBuildersMap; ///< database of modules handled by this instance
  const std::string fType; ///< modules created by this factory
};

#endif
