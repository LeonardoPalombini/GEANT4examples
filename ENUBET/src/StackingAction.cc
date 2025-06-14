#include "StackingAction.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"

#include "RunThread.hh"
#include "G4RunManager.hh"

#include "G4DynamicParticle.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4VCrossSectionDataSet.hh"

StackingAction::StackingAction() {;}


StackingAction::~StackingAction() {;}


G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track* aTrack ) {
  
  G4ClassificationOfNewTrack result( fUrgent );
  
  // G4ParticleDefinition* particleType = aTrack->GetDefinition();
  
  // RunThread* runThread = static_cast<RunThread*>
  //     (G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  
  // if (particleType == G4OpticalPhoton::OpticalPhoton()) { 
  //   G4double energy=aTrack->GetKineticEnergy();
  //      G4int copyNo = aTrack->GetVolume()->GetCopyNo();
  //   if(aTrack->GetVolume()->GetName()=="Scintillator")
  //     runThread->photonEnergyScinti(energy/eV);
  //   else
  //     runThread->photonEnergyWLS(energy/eV);
  // }

  return result;
}


