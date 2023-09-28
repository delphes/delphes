#include "VState.h"

//
// Final state configuration description
//
// Constructor
VState::VState(GenParticle* gMin, std::vector<Track*>& tLink, std::vector<VState*>& vLink)
{
	fMin = gMin;
	fMotherPID = fMin->PID;;
	ftLink = tLink;
	fvLink = vLink;
	//
	// associate gen particles to tracks
	//
	for (UInt_t i = 0; i < ftLink.size(); i++) {
		GenParticle* gen = (GenParticle*) ftLink[i]->Particle.GetObject();
		fpLink.push_back(gen);
	}
}

//
// Destructor
VState::~VState() {
	ftLink.clear();
	fpLink.clear();
	fvLink.clear();
}

//
// Print 
void VState::Print() {
	//
	// print reaction
	TString name = GetPDGname(fMotherPID);
	std::cout << name << " --> ";
	for (UInt_t i = 0; i < fpLink.size(); i++) std::cout << GetPDGname(fpLink[i]->PID) << " ";
	for (UInt_t i = 0; i < fvLink.size(); i++) std::cout << GetPDGname(fvLink[i]->GetMotherPID()) << " ";
	std::cout << endl;
	//
	// Vertex origin and momenta
	std::cout << name << " (Xv= " << fMin->X << ", Yv= " << fMin->Y << ", Zv= " << fMin->Z << ") --- "
		<< "(Px= " << fMin->Px << ", Py= " << fMin->Py << ", Pz= " << fMin->Pz << ")" << std::endl;
	for (UInt_t i = 0; i < fpLink.size(); i++){
		GenParticle* gp = fpLink[i];
		std::cout << GetPDGname(gp->PID) << " (Xv= " << gp->X << ", Yv= " << gp->Y << ", Zv= " << gp->Z << ") --- "
		<< "(Px= " << gp->Px << ", Py= " << gp->Py << ", Pz= " << gp->Pz << ")" << std::endl;
	}
	for (UInt_t i = 0; i < fvLink.size(); i++) {
		VState* vp = fvLink[i];
		GenParticle* gp = vp->GetMother();
		std::cout << GetPDGname(gp->PID) << " (Xv= " << gp->X << ", Yv= " << gp->Y << ", Zv= " << gp->Z << ") --- "
			<< "(Px= " << gp->Px << ", Py= " << gp->Py << ", Pz= " << gp->Pz << ")" << std::endl;
	}
}

//
// Service routines
TString VState::GetPDGname(Int_t pdg)
{
	TString name;
	TParticlePDG* p_pdg = TDatabasePDG::Instance()->GetParticle(pdg);
	if (p_pdg) name = p_pdg->GetName();
	else name = "???";
	//
	return name;
}
