//source for sensitive detector action: hit energy and position

#include "detector.hh"
#include <math.h>



MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name) {
}

MySensitiveDetector::~MySensitiveDetector() {}

G4bool MySensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory* R0hist) {
	/*G4Track* track = aStep->GetTrack();		//retrieve track

	//track->SetTrackStatus(fStopAndKill);


	G4StepPoint* preStepPoint = aStep->GetPreStepPoint();		//in point of photon to det
	G4StepPoint* postStepPoint = aStep->GetPostStepPoint();		//out point of photon from det

	//G4ThreeVector posPhoton = preStepPoint->GetPosition();	//position of photon entrance in a det
	G4double enPhoton = preStepPoint->GetTotalEnergy();


	const G4VTouchable* touchA = aStep->GetPreStepPoint()->GetTouchable();	//get detector of entry
	const G4VTouchable* touchB = aStep->GetPostStepPoint()->GetTouchable();

	G4VPhysicalVolume* physvolA = touchA->GetVolume();
	G4String nameA = physvolA->GetName();		
	G4VPhysicalVolume* physvolB = touchB->GetVolume();
	G4String nameB = physvolB->GetName();

	G4cout << "Names: " << nameA << "  " << nameB << G4endl;

	if (nameA == "physWorld" && nameB == "physScintF") {

		G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();	//# of evt

		auto man = G4AnalysisManager::Instance();			//fill ntuples

		if (nameB == "physScintR") {
			man->FillNtupleIColumn(0, 0, evt);
			man->FillNtupleDColumn(0, 1, enPhoton);
			man->FillNtupleDColumn(0, 2, enPhoton);
			man->AddNtupleRow(0);
		}

		else if (nameB == "physScintF") {
			man->FillNtupleIColumn(1, 0, evt);
			man->FillNtupleDColumn(1, 1, enPhoton);
			man->FillNtupleDColumn(1, 2, enPhoton);
			man->AddNtupleRow(1);
		}

		else {}
	}
	*/
	
	return true;

}