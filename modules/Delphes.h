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

#ifndef Delphes_h
#define Delphes_h

/** \class Delphes
 *
 *  Main Delphes module.
 *  Controls execution of all other modules.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class DelphesConfReader;
class DelphesReader;

class Delphes: public DelphesModule
{
public:
  explicit Delphes(const char *name = "Delphes");

  void Init() override;
  DelphesFactory *GetFactory() const override;

  void SetConfReader(DelphesConfReader *conf) { fConfReader = conf; }
  void SetReader(DelphesReader *reader);
  DelphesReader *GetReader() const { return fReader; }

  void InitTask();
  void ProcessTask();
  void FinishTask();

  void SetOutputFile(std::string_view outputFile) { fOutputFile = outputFile; }
  const std::string &GetOutputFile() const { return fOutputFile; }

  void Clear();

private:
  std::unique_ptr<DelphesFactory> fDelphesFactory;
  DelphesConfReader *fConfReader{nullptr};
  DelphesReader *fReader{nullptr};

  std::vector<std::pair<std::string, std::unique_ptr<DelphesModule> > > fModules;

  std::string fOutputFile;
};

#endif /* Delphes_h */
