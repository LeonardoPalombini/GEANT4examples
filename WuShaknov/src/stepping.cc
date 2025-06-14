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
	G4VPhysicalVolume* fCollim = detectorConstruction->GetCollimator();
	G4VPhysicalVolume* fMon = detectorConstruction->GetMonitor();
	G4VPhysicalVolume* fScat1 = detectorConstruction->GetScatt1();
	G4VPhysicalVolume* fScat2 = detectorConstruction->GetScatt2();
	//G4VPhysicalVolume* fSrc = detectorConstruction->GetSource();

	G4ParticleTable* table = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = step->GetTrack()->GetDefinition();

	
	//kill neutrinos to clean the simulation

	/*if (particle != G4Gamma::Definition()) {
		G4Track* track = step->GetTrack();
		track->SetTrackStatus(fStopAndKill);
		return;
	}*/


	//G4ThreeVector vertex(0., 0., 0.);
	/*if (volume == fScat2 && particle == G4Gamma::Definition()) {
		G4Track* track = step->GetTrack();
		//track->SetTrackStatus(fStopAndKill);
		if (track->GetTrackID() != 1) vertex = track->GetVertexPosition();
		if (vertex.getZ() > 10. * cm && vertex.getZ() < 15. * cm && sqrt(pow(vertex.getX(), 2.) + pow(vertex.getY(), 2.)) < 0.75 * cm) {
			fEventAction->AddStepNo();
		}
		return;
	}*/

	if (0/*volpost == fScat1*/) {
		G4Track* track = step->GetTrack();
		if (track->GetTotalEnergy() < 0.6 * MeV) {
			G4cout << "scat 1 pol:  " << track->GetPolarization() << G4endl;
		}
	}

	if (0/*volpost == fScat2*/) {
		G4Track* track = step->GetTrack();
		if (track->GetTotalEnergy() < 0.6 * MeV) {
			G4cout << "scat 2 pol:  " << track->GetPolarization() << G4endl;
		}
	}


	//when studying the collimator, use monitor to check angular spread only
	if (/*volpost == fMon*/ 0) {
		G4Track* track = step->GetTrack();
		track->SetTrackStatus(fStopAndKill);
		fEventAction->SetMonPosition(step->GetPostStepPoint()->GetPosition());
		//G4cout << "Pol: " << track->GetPolarization() << G4endl;
		G4EventManager* evMan = G4EventManager::GetEventManager();
		evMan->KeepTheCurrentEvent();
		return;
	}

	//when not studying the collimator, it's assumed to be 100% absorbent to speed up the computation
	/*if (volpost == fCollim) {
		G4Track* track = step->GetTrack();
		track->SetTrackStatus(fStopAndKill);

	}*/


	if (particle == G4Gamma::Definition() && (volume == fScat1 || volume == fCollim)) {
		fEventAction->AddScattStepF();
	}
	else if (particle == G4Gamma::Definition() && (volume == fScat2 || volume == fCollim)) {
		fEventAction->AddScattStepR();
	}
	
	
	if (volume == fDetectors.id[0]) {
		//G4Track* track = step->GetTrack();
		//track->SetTrackStatus(fStopAndKill);
		/*G4EventManager* evMan = G4EventManager::GetEventManager();
		evMan->KeepTheCurrentEvent();*/

		G4double edep = step->GetTotalEnergyDeposit();
		fEventAction->AddEdepR(edep);
		//G4EventManager* evMan = G4EventManager::GetEventManager();
		//if(edep > 150*keV) evMan->KeepTheCurrentEvent();

		//fEventAction->SetStepNo(track->GetCurrentStepNumber());
		//G4cout << "Evt " << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << "   step no " << track->GetCurrentStepNumber() << G4endl;
	}
	else if (volume == fDetectors.id[1]) {

		G4double edep = step->GetTotalEnergyDeposit();
		fEventAction->AddEdepF(edep);
		//G4EventManager* evMan = G4EventManager::GetEventManager();
		//if (edep > 150 * keV) evMan->KeepTheCurrentEvent();

	}
	else return;
	
	
}