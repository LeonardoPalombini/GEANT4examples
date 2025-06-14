#include "MuonSensitiveDetector.hh"
#include "RunThread.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4SystemOfUnits.hh"
/*
//TODO: A.B. not anymore in G4-v11
#if(G4VERSION_NUMBER < 1100)
  #include "g4root.hh"
#endif

//TODO: A.B. include in G4-v11
#if(G4VERSION_NUMBER >= 1100)*/
  #include "G4RootAnalysisManager.hh"

MuonSensitiveDetector::MuonSensitiveDetector(G4String SDname)
  : G4VSensitiveDetector(SDname)
{
  G4cout<<"Creating SD with name: "<<SDname<<G4endl;

}


MuonSensitiveDetector::~MuonSensitiveDetector() {

}

G4bool MuonSensitiveDetector::ProcessHits(G4Step* step,G4TouchableHistory*) {
  RunThread* runThread = static_cast<RunThread*>
    (G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  G4TouchableHandle touchable = step->GetPreStepPoint()->GetTouchableHandle();
  G4int copyNo = touchable->GetCopyNumber(0);
  G4int ReplicaNo = touchable->GetReplicaNumber(0);


G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();

  runThread->MuonEnergyDeposit(copyNo, step->GetTotalEnergyDeposit()/MeV);

  return true;

}
//info_ucm->Print();
//hfile.Write();
