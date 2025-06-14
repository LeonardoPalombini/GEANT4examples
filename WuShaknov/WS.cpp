// *** GEANT4 v11 SOFTWARE - CERN ***
// 
// *** Full simulation of the Wu-Shaknov experiment (Physical Review 77.1 (1950): 136) ***
//
// *** scintillator optical response not simulated ***
// 
//Starting file for the build
//performance: about 35k events/min = 600 events/s

#include <iostream>

#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"

#include "construction.hh"
//#include "physics.hh"
#include "PhysicsList.hh"
#include "action.hh"



int main(int argc, char** argv) {

    G4double th = 0.;
    G4cout << "Insert detector angle (deg):"<<G4endl;
    G4cin >> th;

    auto* rMan = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    rMan->SetUserInitialization(new MyDetConstr(th));
    rMan->SetUserInitialization(new PhysicsList());
    rMan->SetUserInitialization(new MyActInit());
    rMan->Initialize();

    G4cout << "Initialization complete!" << G4endl;

    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();


    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    UImanager->ApplyCommand("/vis/open OGL");
    UImanager->ApplyCommand("/vis/viewer/set/viewpointVector 1 1 1");
    UImanager->ApplyCommand("/vis/drawVolume");
    UImanager->ApplyCommand("/vis/viewer/set/autoRefresh true");
    UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");
    UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
    UImanager->ApplyCommand("/vis/scene/add/axes -0.7 -0.7 -0.7 1 m");
    UImanager->ApplyCommand("/vis/geometry/set/forceAuxEdgeVisible");
    UImanager->ApplyCommand("/vis/modeling/trajectories/create/drawByParticleID");
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set gamma green");
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set e- red");
    UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set e+ blue");


    ui->SessionStart();


    return 0;
}
