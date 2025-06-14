#include "RunAction.hh"
#include "RunThread.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Version.hh"
// TODO: A.B. not anymore in G4-v11
#if (G4VERSION_NUMBER < 1100)
#include "g4root.hh"
#endif
// TODO: A.B. include in G4-v11
#if (G4VERSION_NUMBER >= 1100)
#include "G4RootAnalysisManager.hh"
#endif

RunAction::RunAction() : G4UserRunAction()
{

  G4RunManager::GetRunManager()->SetPrintProgress(1);
  G4RootAnalysisManager *analysisManager = G4RootAnalysisManager::Instance();
}

RunAction::~RunAction()
{
  delete G4RootAnalysisManager::Instance();
}

G4Run *RunAction::GenerateRun()
{
  return (new RunThread);
}

void RunAction::BeginOfRunAction(const G4Run *myRun)
{

  /*  m_ROOT_file = new TFile("ucm.root","RECREATE","ROOT file with histograms");
    if(m_ROOT_file) {
      G4cout << "ROOT file energia ucm is created " << G4endl;
    }
  */
  RunThread *runThread = static_cast<RunThread *>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  G4RootAnalysisManager *analysisManager = G4RootAnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);

  G4String fileName;

  if (myRun->GetRunID() == 0)
  {

    if (!analysisManager->GetFileName().empty())
      fileName = analysisManager->GetFileName();
    else
      fileName = "templateAnalysis";

    analysisManager->CreateNtuple("Demo", "Demo");
    analysisManager->CreateNtupleIColumn("PrimaryPDG");
    analysisManager->CreateNtupleDColumn("PosX");
    analysisManager->CreateNtupleDColumn("PosY");
    analysisManager->CreateNtupleDColumn("PosZ");
    analysisManager->CreateNtupleDColumn("PrimaryEne");
    analysisManager->CreateNtupleDColumn("PrimaryDirX");
    analysisManager->CreateNtupleDColumn("PrimaryDirY");
    analysisManager->CreateNtupleDColumn("PrimaryDirZ");
    //analysisManager->CreateNtupleDColumn("eneAbsTot");
    //analysisManager->CreateNtupleDColumn("eneGapTot");
    // analysisManager->CreateNtupleDColumn("muoEneTot");
    analysisManager->CreateNtupleDColumn("muoEne1");
    analysisManager->CreateNtupleDColumn("muoEne2");
    //analysisManager->CreateNtupleDColumn("EneSpess");

    analysisManager->CreateNtupleIColumn("t0Phi", runThread->GetT0Phi());
    analysisManager->CreateNtupleIColumn("t0Z", runThread->GetT0Z());
    analysisManager->CreateNtupleIColumn("caloR", runThread->GetCaloR());
    analysisManager->CreateNtupleIColumn("caloPhi", runThread->GetCaloPhi());
    analysisManager->CreateNtupleIColumn("caloZ", runThread->GetCaloZ());
    analysisManager->CreateNtupleDColumn("t0UpEDep", runThread->GetT0UpEDep());
    analysisManager->CreateNtupleDColumn("t0DwEDep", runThread->GetT0DwEDep());
    analysisManager->CreateNtupleDColumn("ucmEDep", runThread->GetUcmEDep());

    // analysisManager->CreateNtupleDColumn("eneUCM",runThread->GetUcmEnergyDepositVec());
    //  analysisManager->CreateNtupleDColumn("eneAbs",runThread->GetAbsEnergyDepositVec());
    // analysisManager->CreateNtupleDColumn("eneGap",runThread->GetGapEnergyDepositVec());

    //  analysisManager->CreateNtupleIColumn("totPhoton",runThread->GetTotalPhotonVec());
    //  analysisManager->CreateNtupleIColumn("maxPhoton",runThread->GetMaxPhotonVec());
    analysisManager->FinishNtuple();
  }

  analysisManager->OpenFile(fileName);
}

void RunAction::EndOfRunAction(const G4Run * /*run*/)
{

  G4RootAnalysisManager *analysisManager = G4RootAnalysisManager::Instance();

  analysisManager->Write();
  analysisManager->CloseFile();
}

// Create directories
// analysisManager->SetHistoDirectoryName("histograms");
// analysisManager->SetNtupleDirectoryName("ntuple");
// analysisManager->SetVerboseLevel(4);
// analysisManager->SetFirstHistoId(1);

// Creating histograms
// analysisManager->CreateH1("1","Photon Energy in Scintillator", 100, 350., 650.);
// analysisManager->CreateH1("2","Photon Energy in WLS", 100, 350., 650.);

// analysisManager->CreateH1("1","Total Energy in Scintillator", 1000, 0, 10000);
// analysisManager->CreateH1("2","Total Energy in Iron", 1000, 0, 10000);
// analysisManager->CreateH1("3","Photon Energy in SiPM", 100, 350., 650.);
// analysisManager->CreateH1("4","Total Number of Photons in Sensor Plane", 1000, 0, 100000.);

// analysisManager->CreateH1("5","Number of Photons in SiPM1", 1000, 0, 10000.);
// analysisManager->CreateH1("6","Number of Photons in SiPM2", 1000, 0, 10000.);
// analysisManager->CreateH1("7","Number of Photons in SiPM3", 1000, 0, 10000.);
// analysisManager->CreateH1("8","Number of Photons in SiPM4", 1000, 0, 10000.);
// analysisManager->CreateH1("9","Number of Photons in SiPM5", 1000, 0, 10000.);
// analysisManager->CreateH1("10","Number of Photons in SiPM6", 1000, 0, 10000.);
// analysisManager->CreateH1("11","Number of Photons in SiPM7", 1000, 0, 10000.);
// analysisManager->CreateH1("12","Number of Photons in SiPM8", 1000, 0, 10000.);
// analysisManager->CreateH1("13","Number of Photons in SiPM9", 1000, 0, 10000.);
// analysisManager->CreateH1("14","Number of Photons in SiPM10", 1000, 0, 10000.);
// analysisManager->CreateH1("15","Number of Photons in SiPM11", 1000, 0, 10000.);
// analysisManager->CreateH1("16","Number of Photons in SiPM12", 1000, 0, 10000.);
// analysisManager->CreateH1("17","Number of Photons in SiPM13", 1000, 0, 10000.);
// analysisManager->CreateH1("18","Number of Photons in SiPM14", 1000, 0, 10000.);
// analysisManager->CreateH1("19","Number of Photons in SiPM15", 1000, 0, 10000.);
// analysisManager->CreateH1("20","Number of Photons in SiPM16", 1000, 0, 10000.);
