//GEANT4 simulation of 511keV photons Compton scattering on an aluminum block
//study of scattering probabiility and 90° emission 
//scatterer part of WS-like positronium experiment

#include <iostream>

#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"

#include "construction.hh"
#include "physics.hh"
#include "action.hh"


int main(int argc, char** argv) {

    G4double en = 0.;
    G4cout << "Photon energy (MeV):  " << G4endl;
    G4cin >> en;
    G4int n = 1000;
    G4cout << "Number of photons to shoot:  " << G4endl;
    G4cin >> n;
    G4double t = 3.5;
    G4cout << "Tungsten thickness (mm):  " << G4endl;
    G4cin >> t;

	auto* rMan = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    rMan->SetUserInitialization(new MyDetConstr(t));
    rMan->SetUserInitialization(new MyPhysList());
    rMan->SetUserInitialization(new MyActInit(en));
	rMan->Initialize();
    
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
    UImanager->ApplyCommand("/vis/disable");
    UImanager->ApplyCommand("/tracking/storeTrajectory 0");
    UImanager->ApplyCommand("/run/beamOn "+std::to_string(n));

    ui->SessionStart();


	return 0;
}