#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "G4GDMLParser.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "FTFP_BERT.hh"
#include "FTFP_BERT_HP.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BERT_HP.hh"
#include "G4OpticalPhysics.hh"
// TODO: A.B. optical physics managed differently in G4-v11
#if(G4VERSION_NUMBER >= 1100)
	#include "G4OpticalParameters.hh"
#endif
#include "Randomize.hh"


#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4Version.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " template_multi [-m macro ] [-u UIsession] [-t nThreads]" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  if ( argc > 7 ) {
    PrintUsage();
    return 1;
  }
  
  G4String macro;
  G4String session;
#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
#ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) {
      nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else {
      PrintUsage();
      return 1;
    }
  }  
  
  // Detect interactive mode (if no macro provided) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv,"Win32");
  }

  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  if ( nThreads > 0 ) { 
    runManager->SetNumberOfThreads(nThreads);
  }  
#else
  G4RunManager * runManager = new G4RunManager;
#endif

  //ROOT::EnableThreadSafety();

  // Set mandatory initialization classes
  //
  DetectorConstruction* detConstruction = new DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);

  //G4VModularPhysicsList* physicsList = new FTFP_BERT;
 
  //G4VModularPhysicsList* physicsList = new QGSP_BERT_HP;
G4VModularPhysicsList* physicsList = new FTFP_BERT_HP;
  // TODO: A.B. optical physics managed differently in G4-v11
  #if(G4VERSION_NUMBER >= 1100)
    G4OpticalParameters* opticalPhysics = G4OpticalParameters::Instance();
  #else 
    G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  #endif
  //physicsList->SetDefaultCutValue( 0.00007*mm ); 
  
//physicsList->RegisterPhysics(opticalPhysics);
  // adjust some parameters for the optical physics
  opticalPhysics->SetWLSTimeProfile("exponential");
  // TODO: A.B. optical physics managed differently in G4-v11
  #if(G4VERSION_NUMBER < 1100)
    opticalPhysics->SetScintillationYieldFactor(1.0);
    opticalPhysics->SetScintillationExcitationRatio(0.0);
    opticalPhysics->SetMaxNumPhotonsPerStep(100);
    opticalPhysics->SetMaxBetaChangePerStep(10.0);
    opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);
  #endif
 
  if ( physicsList->GetPhysics( "Optical" ) ) {
    G4cout << "***OPTICAL YES***" << G4endl;
  } else {
    G4cout << "***OPTICAL NO*** " << G4endl;
  }
  
  runManager->SetUserInitialization(physicsList);
    
  ActionInitialization* actionInitialization
     = new ActionInitialization();
  runManager->SetUserInitialization(actionInitialization);
  
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( /*macro.size()*/ 0) {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else  {  
    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }

    UImanager->ApplyCommand("/vis/open OGL");
    UImanager->ApplyCommand("/vis/viewer/set/viewpointVector 1 1 1");
    UImanager->ApplyCommand("/vis/drawVolume");
    UImanager->ApplyCommand("/vis/viewer/set/autoRefresh true");
    UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");
    UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
    UImanager->ApplyCommand("/vis/scene/add/axes 0. 0. 0. 1 m");


    ui->SessionStart();
    delete ui;
  }

  
  {           
      G4VPhysicalVolume* g4wv = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume(); 

      G4GDMLParser parser;
      std::ifstream _ff("geo.gdml");
      if (_ff.is_open()) { 
        _ff.close();
        remove("geo.gdml");
      }
      //parser.Write("geo.gdml",g4wv,true,"http://service-spi.web.cern.ch/service-spi/app/releases/GDML/GDML_3_0_0/schema/gdml.xsd");
      parser.SetOverlapCheck (true);
      parser.Write("geo.gdml",g4wv,true);
    }
  




  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
