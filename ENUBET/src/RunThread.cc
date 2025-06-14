#include "RunThread.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
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

RunThread::RunThread() : G4Run()
{

  const DetectorConstruction *constDetConstr =
      static_cast<const DetectorConstruction *>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  DetectorConstruction *detConstr = const_cast<DetectorConstruction *>(constDetConstr);

  
  //nbOfGaps = 1680; // 420;//detConstr->GetNbOfGaps();
  //nbOfAbsor = detConstr->GetNbOfAbsor();
  //nbOfFibers = detConstr->GetNbOfFibers();
  //absNum = detConstr->GetAbsNum();
  //gapNum = detConstr->GetGapNum();
  //sPlaneNum = detConstr->GetSPlaneNum();
  //sipmNum = detConstr->GetSipmNum();
  //nbOfCalo = (detConstr->GetNbOfXCalo()) * (detConstr->GetNbOfYCalo()) * (detConstr->GetNbOfZCalo()) + (detConstr->GetNbOfXCaloH()) * (detConstr->GetNbOfYCaloH()) * (detConstr->GetNbOfZCaloH());


  N_phi_instr = detConstr->GetNPhiInstr();
  N_z_instr = detConstr->GetNZInstr();

  for (auto phi = 0; phi < N_phi_instr; phi++)
    for (auto z = 0; z < N_z_instr; z++)
      for (auto r = -1; r < 3; r++)
      {
        std::vector<G4int> v_coord = {r, phi, z}; // assign coordinates to t0/UCM

        if (r == -1)
        {
          (t0EnergyDepositMap[0])[v_coord] = 0; // set initial value to current t0 doublet
          (t0EnergyDepositMap[1])[v_coord] = 0;
        }
        else
          ucmEnergyDepositMap[v_coord] = 0; // set initial value to current UCM
      }

  // std::cout << detConstr->GetNbOfXCalo() << std::endl;
  // std::cout << detConstr->GetNbOfYCalo() << std::endl;
  // std::cout << detConstr->GetNbOfZCalo() << std::endl;
  // std::cout << detConstr->GetNbOfXCaloH() << std::endl;
  // std::cout << detConstr->GetNbOfYCaloH() << std::endl;
  // std::cout << detConstr->GetNbOfZCaloH() << std::endl;
}

RunThread::~RunThread()
{
}

/*
void RunThread::PhotonCounter(G4int mother, G4int son)
{

  // G4cout<<" mother "<<mother<<" sPlaneNum "<<sPlaneNum<<" son "<<son<<" sipmNum "<<sipmNum<<G4endl;
  photonCounterArray[mother - sPlaneNum][son - sipmNum]++;
}
*/
/*
void RunThread::GapEnergyDeposit(G4int layer, G4double energy)
{

  gapEnergyDepositVec.at(layer - 1) += energy;
}
*/

void RunThread::MuonEnergyDeposit(G4int copyNo, G4double energy) {
	
	if ( copyNo == 1 ) 
  		fEnergyMu1 += energy;
	else if ( copyNo == 2 ) 
		fEnergyMu2 += energy;
}

/* A.B. redefine this method
void RunThread::UcmEnergyDeposit(G4int ucm,G4double energy) {

  ucmEnergyDepositVec.at(ucm)+=energy;
}
*/

void RunThread::AddUcmEnergyDep(std::vector<G4int> ucm, G4double energy)
{
  ucmEnergyDepositMap[ucm] += energy;
}

void RunThread::Addt0EnergyDep(G4int doub, std::vector<G4int> ucm, G4double energy)
{
  (t0EnergyDepositMap[doub])[ucm] += energy;
}
/*
void RunThread::FacEnergyDeposit(G4int fac, G4double energy)
{

  ucmEnergyDepositVec.at(fac) += energy;
}
*/
/*
void RunThread::AbsEnergyDeposit(G4int layer, G4double energy)
{

  absEnergyDepositVec.at(G4int((layer - absNum) / nbOfAbsor)) += energy;
}
*/

void RunThread::FillPerEvent()
{

  // Read t0 and UCM energy for this event and fill analysis NTuple
  for (auto r = -1; r < 3; r++)
    for (auto phi = 0; phi < N_phi_instr; phi++)
      for (auto z = 0; z < N_z_instr; z++)
      {
        std::vector<G4int> v_coord = {r, phi, z};

        if (r == -1) // t0 doublets
        {
          G4double Et0Up = (t0EnergyDepositMap[0])[v_coord];
          G4double Et0Dw = (t0EnergyDepositMap[1])[v_coord];

          if (Et0Up != 0 || Et0Dw != 0)
          {
            t0Phi.push_back(phi);
            t0Z.push_back(z);

            t0UpEDep.push_back(Et0Up);
            t0DwEDep.push_back(Et0Dw);

            //G4cout << "t0(r,phi,z)[0] = (" << r << ", " << phi << ", " << z << ") - E [MeV] = " << (t0EnergyDepositMap[0])[v_coord] << G4endl;
            //G4cout << "t0(r,phi,z)[1] = (" << r << ", " << phi << ", " << z << ") - E [MeV] = " << (t0EnergyDepositMap[1])[v_coord] << G4endl;
          }
        }
        else // UCM calo
        {
          G4double Eucm = ucmEnergyDepositMap[v_coord];

          if (Eucm != 0)
          {
            caloR.push_back(r);
            caloPhi.push_back(phi);
            caloZ.push_back(z);

            ucmEDep.push_back(Eucm);

            //G4cout << "UCM(r,phi,z) = (" << r << ", " << phi << ", " << z << ") - E [MeV] = " << ucmEnergyDepositMap[v_coord] << G4endl;
          }
        }
      }

/*
  for (G4int i = 0; i < nbOfCalo; i++)
  {
    G4int totOpticalPhoton = 0;

    for (G4int j = 0; j < nbOfFibers; j++)
      totOpticalPhoton += photonCounterArray[i][j];

    totalPhotonVec.push_back(totOpticalPhoton);

    G4int maxPhoton = *std::max_element(photonCounterArray[i], photonCounterArray[i] + nbOfFibers);
    maxPhotonVec.push_back(maxPhoton);
  }
  */

  //  std::cout << fPosX << " " << GetPrimaryPosX() << std::endl;
  G4RootAnalysisManager *analysisManager = G4RootAnalysisManager::Instance();
  analysisManager->FillNtupleIColumn(0, fPrimaryPDG);
  analysisManager->FillNtupleDColumn(1, fPosX);
  analysisManager->FillNtupleDColumn(2, fPosY);
  analysisManager->FillNtupleDColumn(3, fPosZ);
  analysisManager->FillNtupleDColumn(4, fEne);
  analysisManager->FillNtupleDColumn(5, fPrimaryDirX);
  analysisManager->FillNtupleDColumn(6, fPrimaryDirY);
  analysisManager->FillNtupleDColumn(7, fPrimaryDirZ);
  // analysisManager->FillNtupleDColumn(8, fEnergyAbsTot);
  // analysisManager->FillNtupleDColumn(9, fEnergyGapTot); // sta riempiendo con l'energia di tutti i gap sommati
  //  if(fEnergyMu > 2){analysisManager->FillNtupleDColumn(10,fEnergyMu);};
  analysisManager->FillNtupleDColumn(8, fEnergyMu1);
  analysisManager->FillNtupleDColumn(9, fEnergyMu2);
  // analysisManager->FillNtupleDColumn(12, fEneSpessore);

  analysisManager->AddNtupleRow();
  // std::cout<< "Energia tile tot " << fEnergyGapTot << std::endl;
}

void RunThread::Reset()
{

  //  std::cout << "Reset" << std::endl;

  // fEnergyGapTot = 0; // Questo è definito in-line in RunThread.hh
  // fEnergyAbsTot = 0;
  //  fEnergyMu = 0;
  fEnergyMu1 = 0;
  fEnergyMu2 = 0;
  // fEneSpessore = 0;

  /*fPrimaryPDG=0;
  fPosX=0;
  fPosY=0;
  fPosZ=0;
  fEne=0;
  fPrimaryDirX=0;
  fPrimaryDirY=0;
  fPrimaryDirZ=0;
  */
  /*
   gapEnergyDepositVec.clear();
   totalPhotonVec.clear();
   maxPhotonVec.clear();
   absEnergyDepositVec.clear();
   ucmEnergyDepositVec.clear();
   */

  /*
    for (G4int i = 0; i < nbOfGaps; i++)
    {

      gapEnergyDepositVec.push_back(0);
      absEnergyDepositVec.push_back(0);

      //  for(G4int j=0;j<nbOfFibers;j++)
      //  photonCounterArray[i][j]=0;
    }

    for (G4int i = 0; i < 336; i++)
    { // for (G4int i = 0; i < 84; i++) {
      ucmEnergyDepositVec.push_back(0);
    }
    */

  t0Phi.clear();
  t0Z.clear();
  caloR.clear();
  caloPhi.clear();
  caloZ.clear();
  t0UpEDep.clear();
  t0DwEDep.clear();
  ucmEDep.clear();

  // reset t0/UCM values
  for (auto phi = 0; phi < N_phi_instr; phi++)
    for (auto z = 0; z < N_z_instr; z++)
      for (auto r = -1; r < 3; r++)
      {
        std::vector<G4int> v_coord = {r, phi, z}; // assign coordinates to t0/UCM

        if (r == -1)
        {
          (t0EnergyDepositMap[0])[v_coord] = 0;
          (t0EnergyDepositMap[1])[v_coord] = 0;
        }
        else
          ucmEnergyDepositMap[v_coord] = 0; // set initial value to current UCM
      }
}

// void RunThread::PhotonEnergyScinti(G4double energy) {
//  G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();
//  analysisManager->FillH1(1,1.23984193*1000/energy);
//}
// void RunThread::PhotonEnergy(G4double energy) {
// G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();
// analysisManager->FillH1(3,1.23984193*1000/energy);
//}

// G4RunManager* runManager =G4RunManager::GetRunManager();

// const DetectorConstruction* detConstr =
//   dynamic_cast<const DetectorConstruction*>
//   (runManager->GetUserDetectorConstruction());

// G4int totOpticalPhoton=0;

///////////////FACTOR//////////////////////////////
// photonCounterArrayMerged[0]=photonCounterArray[0][26]+photonCounterArray[0][27]+
//     photonCounterArray[0][38]+photonCounterArray[0][39];
// photonCounterArrayMerged[1]=photonCounterArray[0][28]+photonCounterArray[0][29]+
//     photonCounterArray[0][40]+photonCounterArray[0][41];
// photonCounterArrayMerged[2]=photonCounterArray[0][30]+photonCounterArray[0][31]+
//     photonCounterArray[0][42]+photonCounterArray[0][43];
// photonCounterArrayMerged[3]=photonCounterArray[0][32]+photonCounterArray[0][33]+
//     photonCounterArray[0][44]+photonCounterArray[0][45];

// photonCounterArrayMerged[4]=photonCounterArray[0][50]+photonCounterArray[0][51]+
//     photonCounterArray[0][62]+photonCounterArray[0][63];
// photonCounterArrayMerged[5]=photonCounterArray[0][52]+photonCounterArray[0][53]+
//     photonCounterArray[0][64]+photonCounterArray[0][65];
// photonCounterArrayMerged[6]=photonCounterArray[0][54]+photonCounterArray[0][55]+
//     photonCounterArray[0][66]+photonCounterArray[0][67];
// photonCounterArrayMerged[7]=photonCounterArray[0][56]+photonCounterArray[0][57]+
//     photonCounterArray[0][68]+photonCounterArray[0][69];

// photonCounterArrayMerged[8]=photonCounterArray[0][74]+photonCounterArray[0][75]+
//     photonCounterArray[0][86]+photonCounterArray[0][87];
// photonCounterArrayMerged[9]=photonCounterArray[0][76]+photonCounterArray[0][77]+
//     photonCounterArray[0][88]+photonCounterArray[0][89];
// photonCounterArrayMerged[10]=photonCounterArray[0][78]+photonCounterArray[0][79]+
//     photonCounterArray[0][90]+photonCounterArray[0][91];
// photonCounterArrayMerged[11]=photonCounterArray[0][80]+photonCounterArray[0][81]+
//     photonCounterArray[0][92]+photonCounterArray[0][93];

// photonCounterArrayMerged[12]=photonCounterArray[0][98]+photonCounterArray[0][99]+
//     photonCounterArray[0][110]+photonCounterArray[0][111];
// photonCounterArrayMerged[13]=photonCounterArray[0][100]+photonCounterArray[0][101]+
//     photonCounterArray[0][112]+photonCounterArray[0][113];
// photonCounterArrayMerged[14]=photonCounterArray[0][102]+photonCounterArray[0][103]+
//     photonCounterArray[0][114]+photonCounterArray[0][115];
// photonCounterArrayMerged[15]=photonCounterArray[0][104]+photonCounterArray[0][105]+
//     photonCounterArray[0][116]+photonCounterArray[0][117];

// for(G4int i=0;i<16;i++)
//   analysisManager->FillH1(i+5,photonCounterArrayMerged[i]);
// for(G4int i=0;i<144;i++)
//   totOpticalPhoton+=photonCounterArray[0][i];

// for(G4int i=0;i<9;i++) {
//   analysisManager->FillH1(i+5,photonCounterArray[0][i]);
//   totOpticalPhoton+=photonCounterArray[0][i];
// }
// analysisManager->FillH1(4,totOpticalPhoton);
// analysisManager->FillH1(1,gapEnergyDeposit[0]);
// analysisManager->FillH1(2,absEnergyDeposit[0]);

// for(G4int j=0;j<maxCalo;j++) {
//   for(G4int i=0;i<maxFiber;i++)
//     photonCounterArray[j][i]=0;
// }

// for(G4int i=0;i<maxCalo;i++)
//   photonCounterArrayMerged[i]=0;

// for(G4int i=0;i<maxCalo;i++) {
//   gapEnergyDeposit[i]=0;
//   absEnergyDeposit[i]=0;
// }
