#ifndef SortableObject_h
#define SortableObject_h

/** \class SortableObject
 *
 *
 *  $Date: 2008-06-04 13:57:26 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "TRef.h"
#include "TObject.h"
#include "TRefArray.h"

#include "TMath.h"

//---------------------------------------------------------------------------

class CompBase
{
public:
  virtual ~CompBase() { }
  virtual Bool_t IsSortable(const TObject *obj) const { return kTRUE; }
  virtual Int_t Compare(const TObject *obj1, const TObject *obj2) const = 0;
};

//---------------------------------------------------------------------------

class SortableObject: public TObject
{
public:

  Bool_t IsSortable() const { return GetCompare() ? GetCompare()->IsSortable(this) : kFALSE; }
  Int_t Compare(const TObject *obj) const { return GetCompare()->Compare(this, obj); }

  virtual const CompBase *GetCompare() const = 0;

  ClassDef(SortableObject, 1)
};

//---------------------------------------------------------------------------
// Standard Comparison Criteria: E, ET, PT, DeltaR
//---------------------------------------------------------------------------

template <typename T>
class CompE: public CompBase
{
  CompE() {}
public:
  static CompE *Instance()
  {
    static CompE single;
    return &single;
  }

  Int_t Compare(const TObject *obj1, const TObject *obj2) const
  {
    const T *t1 = static_cast<const T*>(obj1);
    const T *t2 = static_cast<const T*>(obj2);
    if(t1->E > t2->E)
      return -1;
    else if(t1->E < t2->E)
      return 1;
    else
      return 0;
  }
};

//---------------------------------------------------------------------------

template <typename T>
class CompPT: public CompBase
{
  CompPT() {}
public:
  static CompPT *Instance()
  {
    static CompPT single;
    return &single;
  }

  Int_t Compare(const TObject *obj1, const TObject *obj2) const
  {
    const T *t1 = static_cast<const T*>(obj1);
    const T *t2 = static_cast<const T*>(obj2);
    if(t1->PT > t2->PT)
      return -1;
    else if(t1->PT < t2->PT)
      return 1;
    else
      return 0;
  }
};

//---------------------------------------------------------------------------

template <typename T>
class CompMomentumPt: public CompBase
{
  CompMomentumPt() {}
public:
  static CompMomentumPt *Instance()
  {
    static CompMomentumPt single;
    return &single;
  }

  Int_t Compare(const TObject *obj1, const TObject *obj2) const
  {
    const T *t1 = static_cast<const T*>(obj1);
    const T *t2 = static_cast<const T*>(obj2);
    if(t1->Momentum.Pt() > t2->Momentum.Pt())
      return -1;
    else if(t1->Momentum.Pt() < t2->Momentum.Pt())
      return 1;
    else
      return 0;
  }
};

//---------------------------------------------------------------------------

template <typename T>
class CompET: public CompBase
{
  CompET() {}
public:
  static CompET *Instance()
  {
    static CompET single;
    return &single;
  }

  Int_t Compare(const TObject *obj1, const TObject *obj2) const
  {
    const T *t1 = static_cast<const T*>(obj1);
    const T *t2 = static_cast<const T*>(obj2);
    if(t1->ET > t2->ET)
      return -1;
    else if(t1->ET < t2->ET)
      return 1;
    else
      return 0;
  }
};

//---------------------------------------------------------------------------

template <typename T1, typename T2>
class CompDeltaR: public CompBase
{
  CompDeltaR(const T2 *obj = 0) : fObj(obj) {}

  Double_t DeltaPhi(Double_t phi1, Double_t phi2)
  {
    Double_t phi = TMath::Abs(phi1 - phi2);
    return (phi <= TMath::Pi()) ? phi : (2.0*TMath::Pi()) - phi;
  }

  Double_t Sqr(Double_t x) { return x*x; }

  Double_t SumSqr(Double_t a, Double_t b)
  {
    Double_t aAbs = TMath::Abs(a);
    Double_t bAbs = TMath::Abs(b);
    if(aAbs > bAbs) return aAbs * TMath::Sqrt(1.0 + Sqr(bAbs / aAbs));
    else return (bAbs == 0) ? 0.0 : bAbs * TMath::Sqrt(1.0 + Sqr(aAbs / bAbs));
  };

  const T2 *fObj;

public:
    static CompDeltaR *Instance(const T2 *obj = 0)
  {
      static CompDeltaR single(obj);
      return &single;
  }

  void SetObject(const T2 *obj) { fObj = obj; }

  Int_t Compare(const TObject *obj1, const TObject *obj2) const
  {
    Double_t eta[3], phi[3], deltaR[2];
    const T1 *t1 = static_cast<const T1*>(obj1);
    const T1 *t2 = static_cast<const T1*>(obj2);

    eta[0] = fObj->Eta;
    phi[0] = fObj->Phi;

    eta[1] = t1->Eta;
    phi[1] = t1->Phi;

    eta[2] = t2->Eta;
    phi[2] = t2->Phi;

    deltaR[0] = SumSqr(TMath::Abs(eta[0] - eta[1]), DeltaPhi(phi[0], phi[1]));
    deltaR[1] = SumSqr(TMath::Abs(eta[0] - eta[2]), DeltaPhi(phi[0], phi[2]));

    if(deltaR[0] < deltaR[1])
      return -1;
    else if(deltaR[0] > deltaR[1])
      return 1;
    else
      return 0;
  }
};

#endif // SortableObject_h


