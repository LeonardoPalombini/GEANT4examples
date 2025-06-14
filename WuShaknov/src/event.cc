//source for single-event routine -> here energy deposition in the scintillators

#include "event.hh"

MyEventAction::MyEventAction(MyRunAction*) {
	fEdepR = 0.;
	fEdepF = 0.;
}

MyEventAction::~MyEventAction() {}

void MyEventAction::BeginOfEventAction(const G4Event*) {
	fEdepR = 0.;
	fEdepF = 0.;
	fmonPosition = G4ThreeVector(0., 0., 0.);
	fStepNo = 0.;
	fScatStepsR = 0;
	fScatStepsF = 0;
}

void MyEventAction::EndOfEventAction(const G4Event*) {
	//G4cout << "Energy dep (R,F) : ( " << fEdepR << " , " << fEdepF << " )" << G4endl;
	G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

	G4AnalysisManager* man = G4AnalysisManager::Instance();

	//MyRunAction* runaction = const_cast<MyRunAction*> (static_cast<const MyRunAction*> (G4RunManager::GetRunManager()->GetUserRunAction()));
	
	if (fEdepR > 240.*keV && fEdepR < 300. * keV && fEdepF > 240.*keV && fEdepF < 300. * keV){
		man->FillNtupleIColumn(0, 0, evt);
		man->FillNtupleDColumn(0, 1, 0.);
		man->FillNtupleDColumn(0, 2, fEdepR);
		man->FillNtupleDColumn(0, 3, fScatStepsR);
		man->AddNtupleRow(0);

		man->FillNtupleIColumn(1, 0, evt);
		man->FillNtupleDColumn(1, 1, 0.);
		man->FillNtupleDColumn(1, 2, fEdepF);
		man->FillNtupleDColumn(1, 3, fScatStepsF);
		man->AddNtupleRow(1);

		G4EventManager* evMan = G4EventManager::GetEventManager();
		evMan->KeepTheCurrentEvent();
	}

	//only if studying the collimator
	if (0) {/* fmonPosition.getZ()>100.4 && fmonPosition.getZ() < 104.6*/
		man->FillNtupleDColumn(2, 0, fmonPosition.getX());
		man->FillNtupleDColumn(2, 1, fmonPosition.getY());
		man->FillNtupleDColumn(2, 2, fmonPosition.getZ());
		man->FillNtupleDColumn(2, 3, atan(sqrt(pow(fmonPosition.getX(),2)+ pow(fmonPosition.getY(), 2))/100.));
		man->AddNtupleRow(2);
		
	}


	/*if (fStepNo != 4) {
		
		G4cout << "det hit" << G4endl;
	}*/

	
}

