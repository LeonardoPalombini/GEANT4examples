#ifndef RunThread_h
#define RunThread_h 1

#include "G4Run.hh"
#include "globals.hh"

// max counter 3x3x3calo matrix with 15x15fibers
const G4int maxCalo = 182;
const G4int maxFiber = 9;

class RunThread : public G4Run
{

public:
  RunThread();
  virtual ~RunThread();

  // void PhotonEnergyScinti(G4double);
  // void PhotonEnergy(G4double);
  void FillPerEvent();
  void Reset();
   void PhotonCounter(G4int, G4int);
  // void GapEnergyDeposit(G4int, G4double);

  /* A.B. redefinition of this method
  void UcmEnergyDeposit(G4int,G4double);
<<<<<<< .mine
  */
  /*Add energy deposit for UCM identified by its coordinates ucm = {r, phi, z}*/
  void AddUcmEnergyDep(std::vector<G4int> ucm, G4double energy);
  /*Add energy deposit for t0 identified by the coordinates of the ucm = {-1, phi, z}
  in which is installed. Energy is saved for each doublet doub*/
  void Addt0EnergyDep(G4int doub, std::vector<G4int> ucm, G4double energy);

  // void FacEnergyDeposit(G4int, G4double);
  //void AbsEnergyDeposit(G4int, G4double);
  //   void MuEnergyDeposit(G4double);
  // void GapEnergyDepositTot(G4double ene) { fEnergyGapTot += ene; };
  //void AbsEnergyDepositTot(G4double ene) { fEnergyAbsTot += ene; };
  //  void MuonEnergyDepositTot(G4double ene){fEnergyMu+=ene;}; // uncommented by F. Bramati
  //void MuonEnergyDeposit1(G4double ene) { fEnergyMu1 += ene; }; // edited by F. Bramati
  //void MuonEnergyDeposit2(G4double ene) { fEnergyMu2 += ene; }; // edited by F. Bramati
  void MuonEnergyDeposit(G4int copyNo, G4double ene); // edited by F. Bramati
  // void SpessoreEneDeposit(G4double ene) { fEneSpessore += ene; };

  //G4double GetAbsEnergyDepositTot(void) { return fEnergyAbsTot; };

  G4int GetPrimaryPDG() { return fPrimaryPDG; };
  void SetPrimaryPosX(G4double pos) { fPosX = pos; };
  G4double GetPrimaryPosX(void) { return fPosX; };
  G4double GetPrimaryPosY(void) { return fPosY; };
  void SetPrimaryPosY(G4double pos) { fPosY = pos; };
  void SetPrimaryPosZ(G4double pos) { fPosZ = pos; };
  void SetPrimaryDirX(G4double dir) { fPrimaryDirX = dir; };
  void SetPrimaryDirY(G4double dir) { fPrimaryDirY = dir; };
  void SetPrimaryDirZ(G4double dir) { fPrimaryDirZ = dir; };
  void SetPrimaryEne(G4double ene) { fEne = ene; };
  void SetPrimaryPDG(G4int val) { fPrimaryPDG = val; };
  // std::vector<G4double> &GetGapEnergyDepositVec() { return gapEnergyDepositVec; }
  // std::vector<G4double> &GetUcmEnergyDepositVec() { return ucmEnergyDepositVec; }
  // std::vector<G4double> &GetFacEnergyDepositVec() { return facEnergyDepositVec; }
  //std::vector<G4double> &GetAbsEnergyDepositVec() { return absEnergyDepositVec; }
  //   std::vector<G4double>& GetMuEnergyDepositVec()  {return muEnergyDepositVec;}
  // std::vector<G4int> &GetTotalPhotonVec() { return totalPhotonVec; }
  // std::vector<G4int> &GetMaxPhotonVec() { return maxPhotonVec; }

  // used by analysis manager
  std::vector<G4int> &GetT0Phi() { return t0Phi; };
  std::vector<G4int> &GetT0Z() { return t0Z; };
  std::vector<G4int> &GetCaloR() { return caloR; };
  std::vector<G4int> &GetCaloPhi() { return caloPhi; };
  std::vector<G4int> &GetCaloZ() { return caloZ; };
  std::vector<G4double> &GetT0UpEDep() { return t0UpEDep; };
  std::vector<G4double> &GetT0DwEDep() { return t0DwEDep; };
  std::vector<G4double> &GetUcmEDep() { return ucmEDep; };

private:
  // total (phi,z) sectors instrumented: used to initialize maps
  G4int N_phi_instr;
  G4int N_z_instr;
  // Following maps are used to gather energy deposits by particles crossing UCMs and t0
  std::map<std::vector<G4int>, G4double> ucmEnergyDepositMap;
  std::map<std::vector<G4int>, G4double> t0EnergyDepositMap[2];
  // These vectors store energy deposition per event to be saved in root file (through analysis manager)
  std::vector<G4int> t0Phi;
  std::vector<G4int> t0Z;
  std::vector<G4int> caloR;
  std::vector<G4int> caloPhi;
  std::vector<G4int> caloZ;
  std::vector<G4double> t0UpEDep;
  std::vector<G4double> t0DwEDep;
  std::vector<G4double> ucmEDep;

   G4int photonCounterArray[maxCalo][maxFiber];
   G4int photonCounterArrayMerged[maxCalo];
  // G4double gapEnergyDeposit[maxCalo];
  //G4double absEnergyDeposit[maxCalo];
  // std::vector<G4double> gapEnergyDepositVec;
  // std::vector<G4double> ucmEnergyDepositVec;
  // std::vector<G4double> facEnergyDepositVec;
  //std::vector<G4double> absEnergyDepositVec;
  // std::vector<G4int> totalPhotonVec;
  // std::vector<G4int> maxPhotonVec;

  G4int fPrimaryPDG;
  G4double fPosX;
  G4double fPosY;
  G4double fPosZ;
  G4double fEne;
  G4double fPrimaryDirX;
  G4double fPrimaryDirY;
  G4double fPrimaryDirZ;
  G4double fEnergyAbsTot;
  G4double fEnergyGapTot;
  // G4double fEnergyMu;  // uncommented by F. Bramati
  G4double fEnergyMu1;
  G4double fEnergyMu2;
  G4int nbOfAbsor;
  G4int absNum;
  /*
  G4double fEneSpessore;
  G4int nbOfGaps;
  G4int nbOfFibers;
  G4int nbOfCalo;
  G4int gapNum;
  G4int sPlaneNum;
  G4int sipmNum;
  */
};

#endif
