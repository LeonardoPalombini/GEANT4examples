//source for single-step routine -> here deposited energy integration

#include "stepping.hh"


MySteppingAction::MySteppingAction(MyEventAction* eventAction) {
	fEventAction = eventAction;
}

MySteppingAction::~MySteppingAction() {}

void MySteppingAction::UserSteppingAction(const G4Step* step) {

	G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
	G4VPhysicalVolume* volpost = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
	const MyDetConstr* detectorConstruction = static_cast<const MyDetConstr*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
	PhysDets fDetectors = detectorConstruction->GetPhysDetectors();

	G4ParticleTable* table = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = step->GetTrack()->GetDefinition();


	if ( volume == fDetectors.id ) {

		G4double edep = step->GetTotalEnergyDeposit();
		fEventAction->AddEdep(edep);

		return;
	}

	else return;


}