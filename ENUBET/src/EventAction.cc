#include "EventAction.hh"
#include "RunThread.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "G4Version.hh"
//TODO: A.B. not anymore in G4-v11
#if(G4VERSION_NUMBER < 1100)
  #include "g4root.hh"
#endif
//TODO: A.B. include in G4-v11
#if(G4VERSION_NUMBER >= 1100)
  #include "G4RootAnalysisManager.hh"
#endif
#include <iomanip>

EventAction::EventAction() : G4UserEventAction()
{
}

EventAction::~EventAction() {
}

void EventAction::BeginOfEventAction(const G4Event*) {



  RunThread* runThread = static_cast<RunThread*>
    (G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  runThread->Reset();

}

void EventAction::EndOfEventAction(const G4Event* event) {

  RunThread* runThread = static_cast<RunThread*>
  (G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  runThread->FillPerEvent();

  G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();


  G4int eventID = event->GetEventID();
  G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    //G4cout << "---> End of event: " << eventID << G4endl;

  }
  //if ( rootSaver )
         //{
                 //Retrieve digits collection
                 //G4int digiCollID = digiManager->GetDigiCollectionID( digitsCollName );
                 //const SiDigiCollection* digits = static_cast<const SiDigiCollection*>( digiManager->GetDigiCollection(digiCollID) );
                 //Retrieve hits collections
                 //G4HCofThisEvent* hitsCollections = anEvent->GetHCofThisEvent();
                 //SiHitCollection* hits = 0;
                 //if ( hitsCollections )
                 //{
                         //hits = static_cast<SiHitCollection*>( hitsCollections->GetHC(hitsCollID) );
                 //}
                 //Get Postion and Momentum of primary
                 //This is needed to store in ntuple info @ z=0
                 //const G4ThreeVector& pos = anEvent->GetPrimaryVertex()->GetPosition();
                 //const G4ThreeVector& mom = anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
                 //rootSaver->AddEvent(hits,digits,pos,mom);
      //}
}
