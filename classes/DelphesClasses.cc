
/**
 *
 *  Definition of classes to be stored in the root tree.
 *  Function CompareXYZ sorts objects by the variable XYZ that MUST be
 *  present in the data members of the root tree class of the branch.
 *
 *  $Date: 2008-06-04 13:57:24 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"

#include "classes/DelphesFactory.h"
#include "classes/SortableObject.h"

CompBase *GenParticle::fgCompare = 0;
CompBase *Photon::fgCompare = CompPT<Photon>::Instance();
CompBase *Electron::fgCompare = CompPT<Electron>::Instance();
CompBase *Muon::fgCompare = CompPT<Muon>::Instance();
CompBase *Jet::fgCompare = CompPT<Jet>::Instance();
CompBase *Track::fgCompare = CompPT<Track>::Instance();
CompBase *Tower::fgCompare = CompE<Tower>::Instance();
CompBase *Candidate::fgCompare = CompMomentumPt<Candidate>::Instance();

//------------------------------------------------------------------------------

TLorentzVector GenParticle::P4()
{
  TLorentzVector vec;
  vec.SetPxPyPzE(Px, Py, Pz, E);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector MissingET::P4()
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(MET, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector Photon::P4()
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector Electron::P4()
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector Muon::P4()
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector Jet::P4()
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, Mass);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector Track::P4()
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector Tower::P4()
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(ET, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

Candidate::Candidate() :
  PID(0), Status(0), M1(-1), M2(-1), D1(-1), D2(-1),
  Charge(0), Mass(0.0),
  IsPU(0), IsConstituent(0),
  BTag(0), TauTag(0), Eem(0.0), Ehad(0.0),
  DeltaEta(0.0), DeltaPhi(0.0),
  Momentum(0.0, 0.0, 0.0, 0.0),
  Position(0.0, 0.0, 0.0, 0.0),
  Area(0.0, 0.0, 0.0, 0.0),
  fFactory(0),
  fArray(0)
{
  Edges[0] = 0.0;
  Edges[1] = 0.0;
  Edges[2] = 0.0;
  Edges[3] = 0.0;
}

//------------------------------------------------------------------------------

void Candidate::AddCandidate(Candidate *object)
{
  if(!fArray) fArray = fFactory->NewArray();
  fArray->Add(object);
}

//------------------------------------------------------------------------------

TObjArray *Candidate::GetCandidates()
{
  if(!fArray) fArray = fFactory->NewArray();
  return fArray;
}

//------------------------------------------------------------------------------

Bool_t Candidate::Overlaps(const Candidate *object) const
{
  const Candidate *candidate;

  if(object->GetUniqueID() == GetUniqueID()) return kTRUE;

  if(fArray)
  {
    TIter it(fArray);
    while((candidate = static_cast<Candidate *>(it.Next())))
    {
      if(candidate->Overlaps(object)) return kTRUE;
    }
  }

  if(object->fArray)
  {
    TIter it(object->fArray);
    while((candidate = static_cast<Candidate *>(it.Next())))
    {
      if(candidate->Overlaps(this)) return kTRUE;
    }
  }

  return kFALSE;
}


//------------------------------------------------------------------------------

TObject *Candidate::Clone(const char *newname) const
{
  Candidate *object = fFactory->NewCandidate();
  Copy(*object);
  return object;
}

//------------------------------------------------------------------------------

void Candidate::Copy(TObject &obj) const
{
  Candidate &object = static_cast<Candidate &>(obj);
  Candidate *candidate;

  object.PID = PID;
  object.Status = Status;
  object.M1 = M1;
  object.M2 = M2;
  object.D1 = D1;
  object.D2 = D2;
  object.Charge = Charge;
  object.Mass = Mass;
  object.IsPU = IsPU;
  object.IsConstituent = IsConstituent;
  object.BTag = BTag;
  object.TauTag = TauTag;
  object.Eem = Eem;
  object.Ehad = Ehad;
  object.Edges[0] = Edges[0];
  object.Edges[1] = Edges[1];
  object.Edges[2] = Edges[2];
  object.Edges[3] = Edges[3];
  object.DeltaEta = DeltaEta;
  object.DeltaPhi = DeltaPhi;
  object.Momentum = Momentum;
  object.Position = Position;
  object.Area = Area;

  object.fFactory = fFactory;
  object.fArray = 0;

  if(fArray && fArray->GetEntriesFast() > 0)
  {
    TIter itArray(fArray);
    TObjArray *array = object.GetCandidates();
    while((candidate = static_cast<Candidate *>(itArray.Next())))
    {
      array->Add(candidate);
    }
  }
}

//------------------------------------------------------------------------------

void Candidate::Clear(Option_t* option)
{
  SetUniqueID(0);
  ResetBit(kIsReferenced);
  PID = 0;
  Status = 0;
  M1 = -1; M2 = -1; D1 = -1; D2 = -1;
  Charge = 0;
  Mass = 0.0;
  IsPU = 0;
  IsConstituent = 0;
  BTag = 0;
  TauTag = 0;
  Eem = 0.0;
  Ehad = 0.0;
  Edges[0] = 0.0;
  Edges[1] = 0.0;
  Edges[2] = 0.0;
  Edges[3] = 0.0;
  DeltaEta = 0.0;
  DeltaPhi = 0.0;
  Momentum.SetXYZT(0.0, 0.0, 0.0, 0.0);
  Position.SetXYZT(0.0, 0.0, 0.0, 0.0);
  Area.SetXYZT(0.0, 0.0, 0.0, 0.0);
  fArray = 0;
}
