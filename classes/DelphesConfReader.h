#ifndef classes_DelphesConfReader_h
#define classes_DelphesConfReader_h

/** \class DelphesConfReader
 *
 *  Pure virtual base class handling input card parsing
 *
 *  \author L. Forthomme - AGH, Krakow
 *
 */

#include <string_view>

class DelphesParameters;

class DelphesConfReader
{
public:
  DelphesConfReader() = default;

  virtual void ReadFile(std::string_view fileName) = 0;
  virtual const DelphesParameters &Parameters() const = 0;
};

//------------------------------------------------------------------------------

#endif
