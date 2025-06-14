#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TTree.h"

class G4VPrimaryGenerator;
class G4ParticleGun;
class G4Event;
class PrimaryGeneratorMessenger;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();
  virtual ~PrimaryGeneratorAction();
  virtual void GeneratePrimaries(G4Event *event);

  void SetProduction(G4String val) { fProduction = val; }
  void SetParticle(G4String val) { fParticle = val; }
  void SetEnergy(G4double val) { fEnergy = val; }
  void SetTiltAngle(G4double val) { fTiltAngle = val; }
  void SetRotCenter(G4double val) {fRotCenter = val; }
  void SetBeamType(G4String val) { fBeamType = val; }
  void SetBeamSeed(G4int val) { fBeamSeed = val; }
  void SetBeamProfileFile(G4String val) { fBeamProfileFile = val; }

  const std::vector<G4double>& GetMuMomX()const { return bx; }
  const std::vector<G4double>& GetMuMomY()const { return by; }
  const std::vector<G4double>& GetMuMomZ()const { return bz; }

private:
  //G4VPrimaryGenerator *InitializeGPS();
  G4VPrimaryGenerator *gun;

  // G4VPrimaryGenerator* fParticleGPS;
  G4ParticleGun *fParticleGun;
  PrimaryGeneratorMessenger *gunMessenger; // messenger of this class
  G4String fProduction;                    // flag for a rndm impact point
  G4String fParticle;                      // the particle to be shot
  G4String fBeamType;                      // flag to specify if random beam or profile from histograms
  G4int fBeamSeed;                         // specify seed to extract beam profile from histos
  G4String fBeamProfileFile;               // root file containing histograms with beam profile info

  G4double fEnergy; // particle energy
	
  G4double fTiltAngle; // beam tilt angle in rad
  G4double fRotCenter; // position of rotation center

  TFile *fBeam;
  TH2D *hBSpot; // beam profile from silicon trackers
  //  TH1D *hBTheta; // beam theta divergence
  //  TH1D *hBPhi;   // beam phi divergence
  TH2D *hBangDiv; // beam theta vs phi divergence


  //my addons for muon import
  G4int evtNo=0;
  TTree* muons=nullptr;

  G4int eve = 0;
  G4double muX=0., muY = 0., muZ = 0., muPX = 0., muPY = 0., muPZ = 0.;

  std::vector<G4double> bx;
  std::vector<G4double> by;
  std::vector<G4double> bz;

};

#endif
