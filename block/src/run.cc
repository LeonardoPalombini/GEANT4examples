//source for the run action: create root file

#include "run.hh"


MyRunAction::MyRunAction() {

	auto man = G4AnalysisManager::Instance();
	man->SetVerboseLevel(1);
	man->SetNtupleMerging(true);		//for MT, merge data structures at the end

	man->CreateNtuple("detHits", "detHits");
	man->CreateNtupleIColumn("fEvent");
	man->CreateNtupleDColumn("fEnergy");
	man->FinishNtuple(0);

}


void MyRunAction::BeginOfRunAction(const G4Run* run) {
	auto man = G4AnalysisManager::Instance();

	man->OpenFile("efficiency.root");

	G4cout << "Using " << man->GetType() << G4endl;
}


void MyRunAction::EndOfRunAction(const G4Run* arun) {
	G4AnalysisManager* man = G4AnalysisManager::Instance();

	man->Write();
	man->CloseFile();

}

MyRunAction::~MyRunAction() {}