#ifndef classes_DelphesTCLConfReader_h
#define classes_DelphesTCLConfReader_h

/** \class DelphesTCLConfReader
 *
 *  Class handling TCL input card parsing
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesParameters.h"

#include <tcl/tcl.h>

//------------------------------------------------------------------------------

class DelphesConfReader
{
public:
  DelphesConfReader() = default;

  virtual void ReadFile(std::string_view fileName) = 0;
  virtual const DelphesParameters &Parameters() const = 0;
};

//------------------------------------------------------------------------------

class DelphesTCLConfReader: public DelphesConfReader
{
public:
  DelphesTCLConfReader();

  void ReadFile(std::string_view fileName) override;
  const DelphesParameters &Parameters() const override { return fParams; }

  void SetModuleType(std::string_view moduleName, std::string_view moduleType);

private:
  void ParseValue(const Tcl_Obj *tclObject,
    std::string_view tclName, std::string_view keyName,
    DelphesParameters &delphesParams) const;
  void ParseParameters(std::string_view keyName, DelphesParameters &parametersBlock) const;
  std::vector<std::string> ChildrenList(std::string_view namespaceName) const;
  std::string Run(std::string_view command) const; //FIXME for debugging only!
  std::string TrimmedName(std::string_view keyName) const;

  struct TclInterpDeleter
  {
    void operator()(Tcl_Interp *tclInterp) const { Tcl_DeleteInterp(tclInterp); }
  };
  const std::unique_ptr<Tcl_Interp, TclInterpDeleter> fTclInterp; //!
  std::unordered_map<std::string, std::string> fModuleTypes;
  DelphesParameters fParams;
};

#endif
