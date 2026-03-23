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

#ifndef classes_DelphesParameters_h
#define classes_DelphesParameters_h

#include <yaml-cpp/node/convert.h>
#include <yaml-cpp/yaml.h>

class DelphesParameters
{
public:
  DelphesParameters();
  explicit DelphesParameters(const YAML::Node &); ///< build a parameters block from a YAML structure
  DelphesParameters(const DelphesParameters &);
  virtual ~DelphesParameters() = default;

  DelphesParameters &operator=(const DelphesParameters &); ///< copy assignment operator
  DelphesParameters &operator+=(const DelphesParameters &); ///< add parameters to this parameter object
  DelphesParameters operator+(const DelphesParameters &) const; ///< parameters concatenation operator

  static DelphesParameters loadYAML(const std::string &); ///< load the configuration from a YAML file

  DelphesParameters &Feed(const std::string &); ///< feed a string to steer this parameter block

  const YAML::Node &node() const { return fNode; } ///< retrieve the base node of this parameters container
  bool operator==(const DelphesParameters &oth) const { return oth.fNode == fNode; } ///< equality operator
  friend std::ostream &operator<<(std::ostream &, const DelphesParameters &); ///< human-readable dump of parameters container

  virtual bool Empty() const; ///< does this parameter block contain any parameters?

  std::vector<std::string> Keys() const; ///< list of all parameter keys handled
  DelphesParameters &Erase(const std::string &key); ///< remove a given parameter according to its key
  std::string GetString(const std::string &key) const; ///< string-value of a parameter according to its key

  /// does this parameters block contain a given key?
  template <typename T>
  bool Has(const std::string &key) const
  {
    if(YAML::Node node = fNode[key]; node.IsDefined())
    {
      try
      {
        node.as<T>();
        return true;
      }
      catch(...)
      {
        return false;
      }
    }
    return false;
  }
  /// list of parameter keys with a given type
  template <typename T>
  std::vector<std::string> KeysOf() const
  {
    std::vector<std::string> keysOf;
    for(const std::string &key : Keys())
      if(Has<T>(key))
        keysOf.emplace_back(key);
    return keysOf;
  }
  /// retrieve a parameter with the given key
  template <typename T>
  T Get(const std::string &key, const T &defaultValue = T()) const
  {
    if(YAML::Node node = fNode[key]; node.IsDefined())
      return node.as<T>(defaultValue);
    return defaultValue;
  }
  /// retrieve a parameter with the given key and cast it to a given type
  /// \tparam T expected a stored parameter type
  /// \tparam U type to cast the parameter to
  template <typename T, typename U>
  U GetAs(const std::string &key, const U &defaultValue = U()) const
  {
    return static_cast<U>(Get<T>(key, static_cast<T>(defaultValue)));
  }
  /// set a parameter value
  template <typename T>
  DelphesParameters &Set(const std::string &key, const T &value)
  {
    fNode[key] = value;
    return *this;
  }

private:
  YAML::Node fNode{YAML::NodeType::Map}; ///< base node for the parameter collection
};
static const DelphesParameters kDefaultParameters = DelphesParameters().Set("invalid", true);

using char_array = const char *;

template <>
inline bool DelphesParameters::Has<char_array>(const std::string &key) const
{
  return Has<std::string>(key);
}

template <>
inline DelphesParameters &DelphesParameters::Set(const std::string &key, const char_array &value)
{
  return Set(key, std::string(value));
}

template <>
bool DelphesParameters::Has<DelphesParameters>(const std::string &key) const;

template <>
bool DelphesParameters::Has<std::vector<DelphesParameters> >(const std::string &key) const;

namespace YAML
{
/// conversion rules between a YAML node and a set of integers
template <>
struct convert<std::set<size_t> >
{
  static Node encode(const std::set<size_t> &);
  static bool decode(const Node &, std::set<size_t> &);
};
/// conversion rules between a YAML node and an unordered map of integers
template <>
struct convert<std::unordered_map<size_t, size_t> >
{
  static Node encode(const std::unordered_map<size_t, size_t> &);
  static bool decode(const Node &, std::unordered_map<size_t, size_t> &);
};
/// conversion rules between a YAML node and an unordered map of strings
template <>
struct convert<std::unordered_map<size_t, std::string> >
{
  static Node encode(const std::unordered_map<size_t, std::string> &);
  static bool decode(const Node &, std::unordered_map<size_t, std::string> &);
};
/// conversion rules between a YAML node and an unordered map of integers
template <>
struct convert<DelphesParameters>
{
  static Node encode(const DelphesParameters &);
  static bool decode(const Node &, DelphesParameters &);
};
} // namespace YAML

#endif // classes_DelphesParameters_h
