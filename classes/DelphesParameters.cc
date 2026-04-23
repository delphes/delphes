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

#include <algorithm>
#include <iostream>

#include "classes/DelphesParameters.h"

std::vector<std::string> splitString(std::string_view str, char delim)
{
  std::vector<std::string> out;
  if(str.empty())
    return out;
  std::istringstream iss(std::string{str});
  size_t counter = 0;
  std::string segment;
  while(std::getline(iss, segment, '\"'))
  {
    if(++counter % 2 == 0)
    {
      if(!segment.empty())
        out.emplace_back(segment);
    }
    else
    {
      std::string token;
      std::istringstream iss2(segment);
      while(std::getline(iss2, token, delim))
        if(!token.empty())
          out.emplace_back(token);
    }
  }
  return out;
}

YAML::Node merge_nodes(const YAML::Node &lhs, const YAML::Node &rhs)
{
  if(!rhs.IsMap()) // empty rhs, return lhs if rhs is invalid, otherwise return lhs
    return rhs.IsNull() ? lhs : rhs;
  if(!lhs.IsMap()) // empty lhs, return rhs
    return rhs;
  if(!lhs.size()) // lhs is a map & rhs is an empty map
    return rhs;
  YAML::Node out = lhs; // same mappings as lhs, merged with rhs
  for(YAML::detail::iterator_value node : out)
  {
    try
    {
      const std::string nodeKey = node.first.as<std::string>();
      if(node.first.IsScalar())
      {
        if(const auto key = node.first.Scalar(); rhs[key])
          out[nodeKey] = merge_nodes(node.second, rhs[key]);
        else
          out[nodeKey] = merge_nodes(node.second, lhs[key]);
        continue;
      }
      out[nodeKey] = node.second;
    }
    catch(const YAML::TypedBadConversion<std::string> &err)
    {
      std::ostringstream message;
      message << "Failed to convert node key into a string. Node key:\n"
              << node.first.as<std::string>() << "\n"
              << "yaml-cpp error: " << err.what();
      throw std::runtime_error(message.str());
    }
  }
  for(const YAML::detail::iterator_value &node : rhs)
  { // add mappings from rhs not already added
    try
    {
      const std::string nodeKey = node.first.as<std::string>();
      if(!node.first.IsScalar())
        if(const auto key = node.first.Scalar(); lhs[node.first])
        {
          out[nodeKey] = merge_nodes(out[key], node.second);
          continue;
        }
      out[nodeKey] = merge_nodes(out[nodeKey], node.second);
    }
    catch(const YAML::TypedBadConversion<std::string> &err)
    {
      std::ostringstream message;
      message << "Failed to convert node key into a string. Node key:\n"
              << node.first.as<std::string>() << "\n"
              << "yaml-cpp error: " << err.what();
      throw std::runtime_error(message.str());
    }
  }
  return out;
}

DelphesParameters::DelphesParameters() = default;

DelphesParameters::DelphesParameters(const YAML::Node &node) : fNode(node) {}

DelphesParameters::DelphesParameters(const DelphesParameters &oth) : fNode(Clone(oth.fNode)) {}

DelphesParameters &DelphesParameters::operator=(const DelphesParameters &oth)
{
  fNode = Clone(oth.fNode);
  return *this;
}

DelphesParameters DelphesParameters::operator+(const DelphesParameters &oth) const
{
  DelphesParameters out = *this;
  out += oth;
  return out;
}

DelphesParameters &DelphesParameters::operator+=(const DelphesParameters &oth)
{
  try
  {
    fNode = merge_nodes(fNode, oth.fNode);
    return *this;
  }
  catch(...)
  {
    throw std::runtime_error("Error encountered while concatenating two parameters list.");
  }
}

DelphesParameters DelphesParameters::LoadYAML(std::string_view input_filename)
{
  return DelphesParameters{YAML::LoadFile(std::string{input_filename})};
}

DelphesParameters &DelphesParameters::Feed(std::string_view key_value)
{
  if(std::vector<std::string> tokens = splitString(key_value, ':'); tokens.size() > 1)
  {
    const auto module_key = tokens.at(0);
    tokens.erase(tokens.begin());
    std::string fullString, separator;
    for(const std::string &token : tokens)
      fullString += separator + token, separator = ":";
    Set(module_key, Get<DelphesParameters>(module_key).Feed(fullString));
  }
  else if(const std::vector<std::string> kv = splitString(key_value, '='); kv.size() == 2)
    Set(kv.at(0), kv.at(1));
  else
  {
    std::ostringstream message;
    message << "Failed to unpack key1[:key2[:...]]=value from input string '" << key_value << "'.";
    throw std::runtime_error(message.str());
  }
  return *this;
}

bool DelphesParameters::Empty() const { return !fNode || fNode.size() == 0; }

DelphesParameters &DelphesParameters::Erase(std::string_view key)
{
  if(fNode[key].IsDefined())
    fNode.remove(key);
  return *this;
}

std::string DelphesParameters::GetString(std::string_view key) const
{
  if(auto node = fNode[key]; node.IsDefined())
  {
    if(node.IsSequence())
    {
      if(Has<std::vector<std::string> >(key))
      {
        std::string values, separator;
        for(const std::string &value : Get<std::vector<std::string> >(key))
          values += separator + value, separator = ", ";
        return values;
      }
      if(Has<std::vector<double> >(key))
      {
        std::string values, separator;
        for(const double &value : Get<std::vector<double> >(key))
          values += separator + std::to_string(value), separator = ", ";
        return values;
      }
      if(Has<std::vector<int> >(key))
      {
        std::string values, separator;
        for(const int &value : Get<std::vector<int> >(key))
          values += separator + std::to_string(value), separator = ", ";
        return values;
      }
    }
    if(Has<int>(key))
      return std::to_string(Get<int>(key));
    if(Has<double>(key))
      return std::to_string(Get<double>(key));
    std::ostringstream os;
    os << node;
    return os.str();
  }
  return "{empty}";
}

std::vector<std::string> DelphesParameters::Keys() const
{
  std::vector<std::string> keys;
  if(Empty())
    return keys;
  try
  {
    std::transform(fNode.begin(), fNode.end(), std::back_inserter(keys), [](const YAML::detail::iterator_value &node) {
      try
      {
        return node.first.template as<std::string>();
      }
      catch(const YAML::TypedBadConversion<std::string> &err)
      {
        std::ostringstream message;
        message << "Failed to convert node key into a string. Node key:\n"
                << node.first << "\n"
                << "yaml-cpp error: " << err.what();
        throw std::runtime_error(message.str());
      }
    });
  }
  catch(const YAML::InvalidNode &exception)
  {
    std::cout << "YAML exception: " << exception.msg << std::endl;
  } // ignore all invalid nodes
  return keys;
}

template <>
bool DelphesParameters::Has<DelphesParameters>(std::string_view key) const
{
  // particular case for parameters, as in yaml-cpp, all parameters are technically "Node"s,
  // hence be built/recast as "DelphesParameters" objects
  if(const YAML::Node node = fNode[key]; node.IsDefined() && node.IsMap())
  {
    if(const DelphesParameters &parameters = node.as<DelphesParameters>(kDefaultParameters);
      parameters == kDefaultParameters)
      return false;
    return true;
  }
  return false;
}

template <>
bool DelphesParameters::Has<std::vector<DelphesParameters> >(std::string_view key) const
{
  // particular case for a vector of parameters, as in yaml-cpp, all parameters are technically "Node"s,
  // hence be built/recast as "DelphesParameters" objects
  if(const YAML::Node node = fNode[key]; node.IsDefined() && (node.IsSequence() || node.IsMap()))
  {
    static const std::vector default_parameters_list(2, kDefaultParameters);
    const std::vector<DelphesParameters> &parameters_collection = node.as<std::vector<DelphesParameters> >(default_parameters_list);
    if(parameters_collection == default_parameters_list)
      return false;
    return std::all_of(parameters_collection.begin(), parameters_collection.end(),
      [](const auto &parameter) { return parameter.node().IsMap(); });
  }
  return false;
}

void DelphesParameters::ThrowInvalidConversion(std::string_view keyName) const
{
  std::ostringstream message;
  message << "Invalid conversion requested for key with name '" << keyName << "'. Full parameters block:\n"
          << *this;
  throw std::runtime_error(message.str());
}

std::ostream &operator<<(std::ostream &os, const DelphesParameters &params) { return os << params.fNode; }

namespace YAML
{
Node convert<std::set<size_t> >::encode(const std::set<size_t> &value)
{
  const std::vector vectorView(value.begin(), value.end());
  return convert<std::vector<size_t> >::encode(vectorView);
}

bool convert<std::set<size_t> >::decode(const Node &node, std::set<size_t> &value)
{
  if(!node)
    return false;
  const std::vector<size_t> vectorView = node.as<std::vector<size_t> >();
  value = std::set(vectorView.begin(), vectorView.end());
  return true;
}

Node convert<std::unordered_map<size_t, size_t> >::encode(const std::unordered_map<size_t, size_t> &value)
{
  const std::vector<std::pair<size_t, size_t> > vectorView(value.begin(), value.end());
  return convert<std::vector<std::pair<size_t, size_t> > >::encode(vectorView);
}

bool convert<std::unordered_map<size_t, size_t> >::decode(const Node &node,
  std::unordered_map<size_t, size_t> &value)
{
  if(!node)
    return false;
  const std::vector<std::pair<size_t, size_t> > vector = node.as<std::vector<std::pair<size_t, size_t> > >();
  value = std::unordered_map(vector.begin(), vector.end());
  return true;
}

Node convert<std::unordered_map<size_t, std::string> >::encode(const std::unordered_map<size_t, std::string> &value)
{
  const std::vector<std::pair<size_t, std::string> > vectorView(value.begin(), value.end());
  return convert<std::vector<std::pair<size_t, std::string> > >::encode(vectorView);
}

bool convert<std::unordered_map<size_t, std::string> >::decode(const Node &node,
  std::unordered_map<size_t, std::string> &value)
{
  if(!node)
    return false;
  const std::vector<std::pair<size_t, std::string> > vector = node.as<std::vector<std::pair<size_t, std::string> > >();
  value = std::unordered_map(vector.begin(), vector.end());
  return true;
}

Node convert<DelphesParameters>::encode(const DelphesParameters &params)
{
  if(params.Empty())
    return Node(NodeType::Map);
  try
  {
    return params.node();
  }
  catch(...)
  {
    throw std::runtime_error("Error encountered while encoding parameters list into a YAML node.");
  }
}

bool convert<DelphesParameters>::decode(const Node &node, DelphesParameters &params)
{
  if(!node)
    return false;
  params = DelphesParameters(node);
  return true;
}
} // namespace YAML
