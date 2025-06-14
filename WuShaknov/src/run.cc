//source for the run action: create root file

#include "run.hh"


MyRunAction::MyRunAction() {
	//G4cout << "Set" << G4endl;
	auto man = G4AnalysisManager::Instance();
	man->SetVerboseLevel(1);
	man->SetNtupleMerging(true);		//for MT, merge data structures at the end

	man->CreateNtuple("HitsR", "HitsR");
	man->CreateNtupleIColumn("fEvent");
	man->CreateNtupleDColumn("fEnergy");
	man->CreateNtupleDColumn("fEdep");
	man->CreateNtupleDColumn("fNScatt");
	man->FinishNtuple(0);

	man->CreateNtuple("HitsF", "HitsF");
	man->CreateNtupleIColumn("fEvent");
	man->CreateNtupleDColumn("fEnergy");
	man->CreateNtupleDColumn("fEdep");
	man->CreateNtupleDColumn("fNScatt");
	man->FinishNtuple(1);

	man->CreateNtuple("Mon","Mon");
	man->CreateNtupleDColumn("fMonX");
	man->CreateNtupleDColumn("fMonY");
	man->CreateNtupleDColumn("fMonZ");
	man->CreateNtupleDColumn("fMonTh");
	man->FinishNtuple(2);

}

MyRunAction::~MyRunAction() {}

void MyRunAction::BeginOfRunAction(const G4Run* run) {
	auto man = G4AnalysisManager::Instance();

	man->OpenFile("RUN.root");

	G4cout << "Using " << man->GetType() << G4endl;

	iFEP = 0;
}


void MyRunAction::EndOfRunAction(const G4Run* arun) {
	G4AnalysisManager* man = G4AnalysisManager::Instance();
	man->Write();
	man->CloseFile();

	//G4cout << "FEP events = " << iFEP << G4endl;

}