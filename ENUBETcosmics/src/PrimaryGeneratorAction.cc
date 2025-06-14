#include "PrimaryGeneratorAction.hh"
//#include "PrimaryGeneratorMessenger.hh"
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
#include "G4AnalysisManager.hh"
#include "G4RootAnalysisManager.hh"


#include "TRandom.h"

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction()
{
    /*auto manIN = G4AnalysisManager::Instance();
    manIN->SetVerboseLevel(1);
    manIN->SetNtupleMerging(true);		//for MT, merge data structures at the end

    manIN->CreateNtuple("muonsIN", "muonsIN");
    manIN->CreateNtupleIColumn("fevt");
    manIN->CreateNtupleDColumn("fX");
    manIN->CreateNtupleDColumn("fY");
    manIN->CreateNtupleDColumn("fZ");
    manIN->CreateNtupleDColumn("fPX");
    manIN->CreateNtupleDColumn("fPY");
    manIN->CreateNtupleDColumn("fPZ");
    manIN->FinishNtuple(10);*/
    fBeam = new TFile("muons.root","READ");
    
    muons = (TTree*)fBeam->Get("fMuons");

    muons->SetBranchAddress("fmuX",&muX);
    muons->SetBranchAddress("fmuY", &muY);
    muons->SetBranchAddress("fmuZ", &muZ);
    muons->SetBranchAddress("fmuPX", &muPX);
    muons->SetBranchAddress("fmuPY", &muPY);
    muons->SetBranchAddress("fmuPZ", &muPZ);

    G4int Nent = muons->GetEntries();

    fParticleGun = new G4ParticleGun(1);

    G4ParticleTable* table = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = table->FindParticle("mu-");
    
    fParticleGun->SetParticleEnergy(5.*GeV);
    fParticleGun->SetParticleDefinition(particle);

    evtNo = 0;
    G4cout << "Gun is ready" << G4endl;
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


/*G4VPrimaryGenerator* PrimaryGeneratorAction::InitializeGPS()
{
    G4GeneralParticleSource* gps = new G4GeneralParticleSource();

    // setup details easier via UI commands see gps.mac

    // particle type
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* pion = particleTable->FindParticle("pi+");
    gps->GetCurrentSource()->SetParticleDefinition(pion);

    // set energy distribution
    G4SPSEneDistribution* eneDist = gps->GetCurrentSource()->GetEneDist();
    eneDist->SetEnergyDisType("Mono"); // or gauss
    eneDist->SetMonoEnergy(2.0 * GeV);

    // set position distribution
    G4SPSPosDistribution* posDist = gps->GetCurrentSource()->GetPosDist();
    posDist->SetPosDisType("Beam"); // or Point,Plane,Volume,Beam
    posDist->SetCentreCoords(G4ThreeVector(-2 * m, 140 * cm, 0));
    posDist->SetBeamSigmaInX(0.1 * mm);
    posDist->SetBeamSigmaInY(0.1 * mm);

    // set angular distribution
    G4SPSAngDistribution* angDist = gps->GetCurrentSource()->GetAngDist();
    angDist->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));
    angDist->SetAngDistType("beam2d");
    angDist->SetBeamSigmaInAngX(0.1 * mrad);
    angDist->SetBeamSigmaInAngY(0.1 * mrad);
    angDist->DefineAngRefAxes("angref1", G4ThreeVector(0., 0., 1));

    return gps;
}*/

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    RunThread* runThread = static_cast<RunThread*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    
    muons->GetEntry(evtNo);
    G4ThreeVector pos(muX*cm, muY*cm, muZ*cm);
    G4ThreeVector mom(muPX, muPY, muPZ);
    //G4ThreeVector mom(muPX, muPY, 0.);

    runThread->SetPrimaryPosX(muX);
    runThread->SetPrimaryPosY(muY);
    runThread->SetPrimaryPosZ(muZ);
    runThread->SetPrimaryDirX(muPX);
    runThread->SetPrimaryDirY(muPY);
    runThread->SetPrimaryDirZ(muPZ);

    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(mom);
    fParticleGun->GeneratePrimaryVertex(anEvent);

    evtNo++;
}
