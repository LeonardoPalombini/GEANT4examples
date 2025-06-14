//source for single-event routine -> here energy deposition in the scintillators

#include "event.hh"
#include "Randomize.hh"

MyEventAction::MyEventAction(MyRunAction*) {
	fEdep = 0.;
}

MyEventAction::~MyEventAction() {}

void MyEventAction::BeginOfEventAction(const G4Event*) {
	fEdep = 0.;
	//G4cout << "New event..." << G4endl;

}

void MyEventAction::EndOfEventAction(const G4Event*) {
	//G4cout << "Energy dep:  " << fEdep << G4endl;
	G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

	G4AnalysisManager* man = G4AnalysisManager::Instance();

	//MyRunAction* runaction = const_cast<MyRunAction*> (static_cast<const MyRunAction*> (G4RunManager::GetRunManager()->GetUserRunAction()));
	
	if (fEdep > 10. * keV) {
		man->FillNtupleIColumn(0, 0, evt);
		man->FillNtupleDColumn(0, 1, fEdep /*+ G4RandGauss::shoot(0., fEdep * 0.01 * sqrt(5.035 / fEdep + 0.165))*/);
		man->AddNtupleRow(0);

		/*G4EventManager* evMan = G4EventManager::GetEventManager();
		evMan->KeepTheCurrentEvent();*/
	}

	return;
}