#include "AbsSensitiveDetector.hh"
#include "RunThread.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

AbsSensitiveDetector::AbsSensitiveDetector(G4String SDname)
  : G4VSensitiveDetector(SDname)
{
  G4cout<<"Creating SD with name: "<<SDname<<G4endl;
}

AbsSensitiveDetector::~AbsSensitiveDetector()  { 
}

G4bool AbsSensitiveDetector::ProcessHits(G4Step* step,G4TouchableHistory*) {  

  RunThread* runThread = static_cast<RunThread*>
    (G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  
  G4TouchableHandle touchable = step->GetPreStepPoint()->GetTouchableHandle();
  G4int copyNo = touchable->GetCopyNumber(0);

  /*
  // if(step->GetTrack()->GetParentID()==0)
  if(step->GetTrack()->GetDefinition()->GetPDGEncoding()==12 ||
     step->GetTrack()->GetDefinition()->GetPDGEncoding()==14 ||
     step->GetTrack()->GetDefinition()->GetPDGEncoding()==16 ||
     step->GetTrack()->GetDefinition()->GetPDGEncoding()==-12 ||
     step->GetTrack()->GetDefinition()->GetPDGEncoding()==-14 ||
     step->GetTrack()->GetDefinition()->GetPDGEncoding()==-16)
    {  
        
      // G4cout<<" ANS POS "<<step->GetPostStepPoint()->GetPosition()/cm<<G4endl;
	G4cout << step->GetTrack()->GetDefinition()->GetParticleName() << " "
	//	     << step->GetTrack()->GetDefinition()->GetPDGEncoding() << " "
	     << step->GetPreStepPoint()->GetMaterial()->GetName()<< " "
	  //		     << step->GetPostStepPoint()->GetMaterial()->GetName()<< " "
	     << "E: "     <<step->GetTrack()->GetKineticEnergy()/eV<< " eV "
	     << "Edep: "   << step->GetTotalEnergyDeposit ()/keV << " keV " 
	       << "Edep:so far "   <<runThread->GetAbsEnergyDepositTot()<< " MeV "
	  //     << "gtime: " <<step->GetTrack()->GetGlobalTime()/ns << " "
	  //    <<  " volume: " <<step->GetTrack()->GetVolume()->GetName()<< " "
	  //     <<  " volume: " <<step->GetTrack()->GetVolume()->GetMotherLogical()->GetName()<< " "
	//     <<  " copy: " <<step->GetTrack()->GetVolume()->GetCopyNo()<< " "
	     <<  " position: " << step->GetPostStepPoint()->GetPosition()/cm << " " 
        ////<<  " step " <<   step->GetStepLength()/um << " "   
	//     <<  " ID: " <<step->GetTrack()->GetTrackID() << " " 
        //<<  " Parent ID: "    <<step->GetTrack()->GetParentID() << " " 
        //// <<  " ltime: " <<step->GetTrack()->GetLocalTime()/ns  << " "
        //// <<  " steps: " <<step->GetTrack()->GetCurrentStepNumber() << " "
	//     << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName ()
	     << G4endl ;

    }
  */

  //runThread->AbsEnergyDeposit(copyNo,step->GetTotalEnergyDeposit()/MeV);
  //runThread->AbsEnergyDepositTot(step->GetTotalEnergyDeposit()/MeV);
  
  return true;
}
