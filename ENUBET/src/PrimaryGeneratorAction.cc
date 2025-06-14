#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GeneralParticleSource.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "RunThread.hh"
#include "G4RunManager.hh"
#include "G4Version.hh"
// TODO: A.B. not anymore in G4-v11
#if (G4VERSION_NUMBER < 1100)
#include "g4root.hh"
#endif
// TODO: A.B. include in G4-v11
#if (G4VERSION_NUMBER >= 1100)
#include "G4RootAnalysisManager.hh"
#endif

#include "TRandom.h"

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction()
{
  // gun = InitializeGPS();
  //  create a messenger for this class
  gunMessenger = new PrimaryGeneratorMessenger(this);

  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  G4ParticleDefinition *particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("pi+");
  fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticlePosition(G4ThreeVector(-50. * cm, 138. * cm, 0. * cm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
  fParticleGun->SetParticleEnergy(1 * GeV);

  // fProduction="no";
  fBeamType = "histo"; // this is the default value, if not specified in the mac file
  fBeamSeed = -1;      // this is the default value, if not specified in the mac file
	
  fTiltAngle = 0.*rad; // this is the default value, if not specified in the mac file
  fRotCenter = 0.*cm; //this is the defaul value, if not specified in the mac file
	
  // get beam profile info
  fBeam = NULL;

  hBSpot = NULL;
  hBangDiv = NULL;
  //  hBTheta = NULL;
  //  hBPhi = NULL;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  // delete fParticleGPS;
  delete gunMessenger;
  // delete gun;

  // delete gunMessenger;
  delete fParticleGun;

  delete fBeam;
  delete hBSpot;
  delete hBangDiv;
  //  delete hBTheta;
  //  delete hBPhi;
}

G4VPrimaryGenerator *PrimaryGeneratorAction::InitializeGPS()
{
  G4GeneralParticleSource *gps = new G4GeneralParticleSource();

  // setup details easier via UI commands see gps.mac

  // particle type
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *pion = particleTable->FindParticle("pi+");
  gps->GetCurrentSource()->SetParticleDefinition(pion);

  // set energy distribution
  G4SPSEneDistribution *eneDist = gps->GetCurrentSource()->GetEneDist();
  eneDist->SetEnergyDisType("Mono"); // or gauss
  eneDist->SetMonoEnergy(2.0 * GeV);

  // set position distribution
  G4SPSPosDistribution *posDist = gps->GetCurrentSource()->GetPosDist();
  posDist->SetPosDisType("Beam"); // or Point,Plane,Volume,Beam
  posDist->SetCentreCoords(G4ThreeVector(-2 * m, 140 * cm, 0));
  posDist->SetBeamSigmaInX(0.1 * mm);
  posDist->SetBeamSigmaInY(0.1 * mm);

  // set angular distribution
  G4SPSAngDistribution *angDist = gps->GetCurrentSource()->GetAngDist();
  angDist->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));
  angDist->SetAngDistType("beam2d");
  angDist->SetBeamSigmaInAngX(0.1 * mrad);
  angDist->SetBeamSigmaInAngY(0.1 * mrad);
  angDist->DefineAngRefAxes("angref1", G4ThreeVector(0., 0., 1));

  return gps;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  // at first event, get the histograms with beam profile info
  if (fBeamType == "histo" && fBeam == NULL)
  {
    if (fBeamSeed != -1)
      gRandom->SetSeed(fBeamSeed);

    fBeam = new TFile(fBeamProfileFile);

    if (fBeam->IsZombie())
    {
      G4Exception("PrimaryGeneratorAction::GeneratePrimaries()", "", FatalException, "Error: root file with beam profile info is bad or not provided!");
    }

    //    hBSpot = (TH2D *)fBeam->Get("hBSpot");
    //    hBTheta = (TH1D *)fBeam->Get("hBTheta");
    //    hBPhi = (TH1D *)fBeam->Get("hBPhi");
    hBSpot = (TH2D *)fBeam->Get("beamSpot");
    hBangDiv = (TH2D *)fBeam->Get("angDiv");
  }

  RunThread *runThread = static_cast<RunThread *>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  /*if(fProduction=="flat")
    {
      G4double x0 = 0.*cm;
      G4double y0 = 0.*cm;
      G4double z0 = 0.*cm;

      //standard simu

      z0 = (G4UniformRand()*7-3.5)*cm;
      y0 = (G4UniformRand()*7-3.5)*cm;
      x0= -7.5*cm;
      stop comment*/

  G4RootAnalysisManager *analysisManager = G4RootAnalysisManager::Instance();

  /*G4GeneralParticleSource *fParticleGPS = (G4GeneralParticleSource *)gun;

  fParticleGPS->GeneratePrimaryVertex(anEvent);

  runThread->SetPrimaryPDG(fParticleGPS->GetParticleDefinition()->GetPDGEncoding());
  runThread->SetPrimaryPosX(fParticleGPS->GetParticlePosition().x());
  runThread->SetPrimaryPosY(fParticleGPS->GetParticlePosition().y());
  runThread->SetPrimaryPosZ(fParticleGPS->GetParticlePosition().z());
  runThread->SetPrimaryDirX(fParticleGPS->GetParticleMomentumDirection().x());
  runThread->SetPrimaryDirY(fParticleGPS->GetParticleMomentumDirection().y());
  runThread->SetPrimaryDirZ(fParticleGPS->GetParticleMomentumDirection().z());
  runThread->SetPrimaryEne(fParticleGPS->GetParticleEnergy());

  // std::cout << fParticleGPS->GetParticlePosition().x() << " " << runThread->GetPrimaryPosX() << std::endl;

  // FIXME: A.B. test
  G4cout << "-------PART GENERATION PDG = " << fParticleGPS->GetParticleDefinition()->GetPDGEncoding() << G4endl;
  G4cout << "ps x = " << fParticleGPS->GetParticlePosition().x() << G4endl;
  G4cout << "ps y = " << fParticleGPS->GetParticlePosition().y() << G4endl;
  G4cout << "ps z = " << fParticleGPS->GetParticlePosition().z() << G4endl;
  G4cout << "dir x = " << fParticleGPS->GetParticleMomentumDirection().x() << G4endl;
  G4cout << "dir y = " << fParticleGPS->GetParticleMomentumDirection().y() << G4endl;
  G4cout << "dir z = " << fParticleGPS->GetParticleMomentumDirection().z() << G4endl;
  G4cout << "E = " << fParticleGPS->GetParticleEnergy() << G4endl;*/

  G4ParticleDefinition *particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(fParticle);
  fParticleGun->SetParticleDefinition(particleDefinition);

  // A.B. Example to generate a flat square beam spot, with beam divergence described by
  // theta angle along y and phi angle along z
  //  G4double dz = 10. * (G4UniformRand() - 0.5);
  //  G4double dy = 10. * (G4UniformRand() - 0.5);

  // S.M. Generate a flat square beam spot with beam divergence described by Theta angle along Y and Phi angle along Z. Their values are sampled from their distribution taken by Silicium Trackers during TB2022
  G4double dz = 0. * cm;
  G4double dy = 0. * cm;
  G4double dx = 0. * cm;
  G4double theta = 0.;
  G4double phi = 0.;
  if (fBeamType == "random")
  {
    dz = 10. * (G4UniformRand() - 0.5);
    dy = 10. * (G4UniformRand() - 0.5);
    theta = 0.5 * (3.14 / 180.) * G4RandGauss::shoot(0, 1);
    phi = 0.5 * (3.14 / 180.) * G4RandGauss::shoot(0, 1);
  }
  else if (fBeamType == "histo")
  {
    hBSpot->GetRandom2(dz, dy);
    hBangDiv->GetRandom2(phi, theta);
  }
	//S.M.: to consider tilted beam angle, it is necessary to rotate both Beam Spot position and Primary Particles directions of fTiltAngle value along z-axis, shifting the coordinates of "fRotCenter" to consider the center of rotation in the origin of World coordinates.
  // Distance between Si1 and front face of the Demonstrator (see LogBook)
  G4double distance = ((4.7 / 2.) + 69.5 + 4.7 + 90.5) * cm;
  dx = -fRotCenter-distance;
	
  //2D Rotation along z-axis with center in the origin of World coordinate system
  G4double dx_rot = 0*cm;
  G4double dy_rot = 0*cm;
	
  dx_rot = dx*cos(fTiltAngle)-(dy*cm)*sin(fTiltAngle);
  dy_rot = dx*sin(fTiltAngle)+(dy*cm)*cos(fTiltAngle);
  //To reset original coordinate system
  dx_rot += fRotCenter;
  dy_rot += 136.5*cm;  //136.5 cm is the height from ground of the border between t0 and t1 layer (considering y-axis). This is the center of the beam spot position
  fParticleGun->SetParticlePosition(G4ThreeVector(dx_rot, dy_rot, dz*cm));
//  G4double norm = sqrt(1 + pow(tan(theta), 2) + pow(tan(phi), 2));
  // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(cos(theta), sin(theta), 0.));
	
  //Same 2D Rotation has to be performed on Primary Particle Directions with angle = - fTiltAngle
  G4double dirX = 1;
  G4double dirY = tan(theta);
  G4double dirZ = tan(phi);
	
  G4double dirX_rot = 0*rad;
  G4double dirY_rot = 0*rad;
  dirX_rot = dirX*cos(fTiltAngle)-dirY*sin(fTiltAngle);
  dirY_rot = dirX*sin(fTiltAngle)+dirY*cos(fTiltAngle);
	
	//DEBUG PART
//	G4cout<<"---------*******--------------******--------------"<<G4endl;
//	G4cout<<"fTiltAngle: "<<fTiltAngle <<G4endl;
//	G4cout<<"dx: "<<dx/cm <<G4endl;
//	G4cout<<"dy: "<<(dy*cm)/cm <<G4endl;
//	G4cout<<"dz: "<<(dz*cm)/cm <<G4endl;
//	G4cout<<"dx_rot: "<<dx_rot/cm <<G4endl;
//	G4cout<<"dy_rot: "<<dy_rot/cm <<G4endl;
//	G4cout<<"dirX: "<<dirX_rot <<G4endl;
//	G4cout<<"dirY: "<<dirY_rot <<G4endl;
//	G4cout<<"dirZ: "<<dirZ <<G4endl;
//	G4cout<<"---------*******--------------******--------------"<<G4endl;
	//
	
  G4double norm = sqrt(pow(dirX_rot,2)+pow(dirY_rot,2)+pow(dirZ,2));
	
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(dirX_rot/norm,dirY_rot/norm,dirZ/norm/*1 / norm, tan(theta) / norm, tan(phi) / norm*/));

  fParticleGun->SetParticleEnergy(fEnergy * GeV);

  fParticleGun->GeneratePrimaryVertex(anEvent);

  runThread->SetPrimaryPDG(fParticleGun->GetParticleDefinition()->GetPDGEncoding());
  runThread->SetPrimaryPosX(fParticleGun->GetParticlePosition().x() / cm);
  runThread->SetPrimaryPosY(fParticleGun->GetParticlePosition().y() / cm);
  runThread->SetPrimaryPosZ(fParticleGun->GetParticlePosition().z() / cm);
  runThread->SetPrimaryDirX(fParticleGun->GetParticleMomentumDirection().x());
  runThread->SetPrimaryDirY(fParticleGun->GetParticleMomentumDirection().y());
  runThread->SetPrimaryDirZ(fParticleGun->GetParticleMomentumDirection().z());
  runThread->SetPrimaryEne(fParticleGun->GetParticleEnergy());

  // FIXME: A.B. test
  /*G4cout << "-------PART GENERATION PDG = " << fParticleGun->GetParticleDefinition()->GetPDGEncoding() << G4endl;
  G4cout << "ps x = " << fParticleGun->GetParticlePosition().x() << G4endl;
  G4cout << "ps y = " << fParticleGun->GetParticlePosition().y() << G4endl;
  G4cout << "ps z = " << fParticleGun->GetParticlePosition().z() << G4endl;
  G4cout << "dir x = " << fParticleGun->GetParticleMomentumDirection().x() << G4endl;
  G4cout << "dir y = " << fParticleGun->GetParticleMomentumDirection().y() << G4endl;
  G4cout << "dir z = " << fParticleGun->GetParticleMomentumDirection().z() << G4endl;
  G4cout << "E = " << fParticleGun->GetParticleEnergy() << G4endl;*/
}
