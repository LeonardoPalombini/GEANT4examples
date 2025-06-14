#include "SipmSensitiveDetector.hh"
#include "RunThread.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4SystemOfUnits.hh"

SipmSensitiveDetector::SipmSensitiveDetector(G4String SDname)
  : G4VSensitiveDetector(SDname)
{
  G4cout<<"Creating SD with name: "<<SDname<<G4endl;
}

SipmSensitiveDetector::~SipmSensitiveDetector() { 
}

G4bool SipmSensitiveDetector::ProcessHits(G4Step* step,G4TouchableHistory*) {  

  G4TouchableHandle touchable = step->GetPreStepPoint()->GetTouchableHandle();
  G4int copyNoSon = touchable->GetCopyNumber(0);
  G4int copyNoMother = touchable->GetCopyNumber(1);
  
  RunThread* runThread = static_cast<RunThread*>
    (G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  /*
 if(step->GetTrack()->GetParentID()==0)
    {    
      
      //  G4cout<<" ANS POS "<<step->GetPostStepPoint()->GetPosition()/cm<<G4endl;
      G4cout << step->GetTrack()->GetDefinition()->GetParticleName() << " "
	     << step->GetTrack()->GetDefinition()->GetPDGEncoding() << " "
	     << step->GetPreStepPoint()->GetMaterial()->GetName()<< " "
	     << step->GetPostStepPoint()->GetMaterial()->GetName()<< " "
	     << "E: "     <<step->GetTrack()->GetKineticEnergy()/eV<< " eV "
	     << "Edep: "   << step->GetTotalEnergyDeposit ()/keV << " keV " 
	     << "gtime: " <<step->GetTrack()->GetGlobalTime()/ns << " "
	     <<  " volume: " <<step->GetTrack()->GetVolume()->GetName()<< " "
	     <<  " volume: " <<step->GetTrack()->GetVolume()->GetMotherLogical()->GetName()<< " "
	     <<  " copy: " <<step->GetTrack()->GetVolume()->GetCopyNo()<< " "
	     <<  " position: " << step->GetPostStepPoint()->GetPosition()/cm << " " 
        //<<  " step " <<   step->GetStepLength()/um << " "   
	     <<  " ID: " <<step->GetTrack()->GetTrackID() << " " 
        <<  " Parent ID: "    <<step->GetTrack()->GetParentID() << " " 
        // <<  " ltime: " <<step->GetTrack()->GetLocalTime()/ns  << " "
        // <<  " steps: " <<step->GetTrack()->GetCurrentStepNumber() << " "
	     << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName ()
	     << G4endl ;

    }*/

/*
  G4ParticleDefinition* particleType = step->GetTrack()->GetDefinition();
  if (particleType == G4OpticalPhoton::OpticalPhoton()) { 
    runThread->PhotonCounter(copyNoMother,copyNoSon);
    // G4double energy=step->GetTrack()->GetKineticEnergy();
    // runThread->PhotonEnergy(energy/eV);
  }
  */
  
  return true;
}
