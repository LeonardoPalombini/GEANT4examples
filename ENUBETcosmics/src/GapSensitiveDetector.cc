#include "GapSensitiveDetector.hh"
#include "RunThread.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "DetectorConstruction.hh"
#include "G4Version.hh"
// TODO: A.B. not anymore in G4-v11
#if (G4VERSION_NUMBER < 1100)
#include "g4root.hh"
#endif
// TODO: A.B. include in G4-v11
#if (G4VERSION_NUMBER >= 1100)
#include "G4RootAnalysisManager.hh"
#endif

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
using namespace std;

GapSensitiveDetector::GapSensitiveDetector(G4String SDname)
    : G4VSensitiveDetector(SDname)
{
  G4cout << "Creating SD with name: " << SDname << G4endl;
}

GapSensitiveDetector::~GapSensitiveDetector()
{
}

void GapSensitiveDetector::Initialize(G4HCofThisEvent *hce)
{
  G4RunManager *runManager = G4RunManager::GetRunManager();

  const DetectorConstruction *detConstr = dynamic_cast<const DetectorConstruction *>(runManager->GetUserDetectorConstruction());

  gapCoord = detConstr->GetGapCoord();
}

G4bool GapSensitiveDetector::ProcessHits(G4Step *step, G4TouchableHistory *)
{
  RunThread *runThread = static_cast<RunThread *>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  G4TouchableHandle touchable = step->GetPreStepPoint()->GetTouchableHandle();
  G4int copyNo = touchable->GetCopyNumber(0);
  // G4int ReplicaNo = touchable->GetReplicaNumber(0);

  G4double energy = step->GetTotalEnergyDeposit() / MeV;

  // G4double tile_ucm[84];
  //  prova per creare ucm con vector di vector
  //  da mettere nel private .hh
  //  std::vector<G4double>& GetucmVec() {return ucmVec;}
  // std::vector<vector<double>> ucmVec;

  // da utilizzare
  // devo fare un vector con un vector come elementi
  // riempio il std::vector<int> v;
  // for (size_t i = 0; i < 84; i++)
  //{
  //  for (size_t j = 6; j < 6; j++)
  //  {
  //    ucmVec[i].push_back(i);
  //  }
  //}
  //
  // G4double a;
  //
  // G4double xposition = step->GetPreStepPoint()->GetPosition().x();
  // G4double yposition = step->GetPreStepPoint()->GetPosition().y();
  // G4double zposition = step->GetPreStepPoint()->GetPosition().z();
  //
  // if (yposition >= -60 && yposition < -30)
  //{
  //  a = -30;
  //}
  // if (yposition >= -30 && yposition < 0)
  //{
  //  a = 0;
  //}
  // if (yposition >= 0 && yposition < 30)
  //{
  //  a = 30;
  //}
  // if (yposition >= 30 && yposition < 60)
  //{
  //  a = 60;
  //  G4cout << a << G4endl;
  //}

  // else a = 0;

  // G4double ypos = 15 - (a - yposition);
  // G4double energydep = (step->GetTotalEnergyDeposit() / MeV) * (0.7547 + pow(ypos, 2) * 0.001537) * 1.5;
  // G4double id = step->GetTrack()->GetTrackID();
  // G4double pdg = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();

  // G4RootAnalysisManager *analysisManager = G4RootAnalysisManager::Instance();
  //  runThread->GapEnergyDeposit(copyNo, step->GetTotalEnergyDeposit() / MeV);
  //  runThread->GapEnergyDepositTot(step->GetTotalEnergyDeposit() / MeV);
  //   runThread->UcmEnergyDeposit((copyNo-1)/5,step->GetTotalEnergyDeposit()/MeV);
  //  runThread->UcmEnergyDeposit((copyNo - 1) / 5, energydep);
  //   runThread->FacEnergyDeposit((copyNo-1)/5,step->GetTotalEnergyDeposit()/MeV);
  //   G4cout << a << "ypos " <<ypos<<" yposition "<<yposition<< G4endl;

  // get tile (gap) coordinates in (r, phi, x)
  std::vector<G4int> vc = gapCoord[copyNo];
  G4int r = vc[0];
  G4int phi = vc[1];
  G4int x = vc[2];
  // we need to transofrm the x layer of the tile into the corrisponding z layer
  // of the UCM (each UCM is made out of 5 tiles in x)
  G4int z = (int)x / 5;
  // if a t0 we need to know if is upward/downward tile
  G4int doublet = x - z * 5; // value either 0 or 1

  std::vector<G4int> ucm = {r, phi, z};

  // Add energy to t0-doublet or UCM calorimeter, depending on the tile radial position
  if (r == -1)
  {
    runThread->Addt0EnergyDep(doublet, ucm, energy);
  }
  else
  {
    runThread->AddUcmEnergyDep(ucm, energy);
  }

  // FIXME: A.B. test
  /*G4cout << "copyNo = " << copyNo << G4endl;
  G4cout << "vc.size() = " << vc.size() << G4endl;
  G4cout << "tile coordinates: R = " << vc[0] << " - Phi = " << vc[1] << " - X = " << vc[2] << G4endl;
  // G4cout << "UCM coordinates: R = " << ucm[0] << " - Phi = " << ucm[1] << " - Z = " << ucm[2] << G4endl;
  G4cout << "name = " << step->GetTrack()->GetVolume()->GetName() << G4endl;
  // G4cout << "pdg = " << pdg << G4endl;
  G4cout << "E [MeV] = " << step->GetTotalEnergyDeposit() / MeV << G4endl;
  // G4cout << "(x,y,z) = (" << xposition << ", " << yposition << ", " << zposition << ")" << G4endl;*/

  return true;
}
// info_ucm->Print();
// hfile.Write();
