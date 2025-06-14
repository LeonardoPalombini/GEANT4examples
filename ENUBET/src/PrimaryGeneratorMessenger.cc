
#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
    PrimaryGeneratorAction *Gun)
    : Action(Gun)
{
  gunDir = new G4UIdirectory("/demo/gun/");
  gunDir->SetGuidance("PrimaryGenerator control");

  ProductionCmd = new G4UIcmdWithAString("/demo/gun/production", this);
  ProductionCmd->SetGuidance("Select type of production");
  ProductionCmd->SetGuidance("  Choice :flat tilt no(default) ");
  ProductionCmd->SetParameterName("production", true);
  ProductionCmd->SetDefaultValue("no");
  ProductionCmd->SetCandidates("no flat tilt");
  ProductionCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  BeamTypeCmd = new G4UIcmdWithAString("/demo/gun/beamType", this);
  BeamTypeCmd->SetGuidance("Select beam type");
  BeamTypeCmd->SetGuidance("  Choice : random histo(default) ");
  BeamTypeCmd->SetParameterName("beamType", true);
  BeamTypeCmd->SetDefaultValue("histo");
  BeamTypeCmd->SetCandidates("random histo");
  BeamTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  ParticleCmd = new G4UIcmdWithAString("/demo/gun/particle", this);
  ParticleCmd->SetGuidance("Select particle (e+, pi+, etc.)");
  ParticleCmd->SetDefaultValue("");

  EnergyCmd = new G4UIcmdWithADouble("/demo/gun/energy", this);
  EnergyCmd->SetGuidance("Select particle energy");
  EnergyCmd->SetDefaultValue(0);

  BeamSeedCmd = new G4UIcmdWithAnInteger("/demo/gun/beamSeed", this);
  BeamSeedCmd->SetGuidance("Set seed for generating beam profile from histos (to be used when providing a beamProfileFile). If 0 each simulation run will have independent events.");
  BeamSeedCmd->SetDefaultValue(0);

  BeamProfileFileCmd = new G4UIcmdWithAString("/demo/gun/beamProfileFile", this);
  BeamProfileFileCmd->SetGuidance("Select root file containing histos describing beam profile (beam spot, divergence, etc.)");
  BeamProfileFileCmd->SetDefaultValue("beamProfile.root");
	
  //S.M. : Tilted beam Command
  TiltedBeamCmd = new G4UIcmdWithADouble("/demo/gun/beamTiltAngle", this);
  TiltedBeamCmd->SetGuidance("Select beam tilt angle (in rad)");
  TiltedBeamCmd->SetDefaultValue(0);
	
  //S.M.: Rotation center distance Command (position of rotation center from Front Face of Demo, X=0)
  RotCenterCmd = new G4UIcmdWithADouble("/demo/gun/rotCenterDist", this);
  RotCenterCmd->SetGuidance("Select distance of rotation center from front face of demonstrator (in mm)");
  RotCenterCmd->SetDefaultValue(0);
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete ProductionCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(
    G4UIcommand *command, G4String newValue)
{
  if (command == ProductionCmd)
  {
    Action->SetProduction(newValue);
  }
  else if (command == BeamTypeCmd)
  {
    Action->SetBeamType(newValue);
  }
  else if (command == ParticleCmd)
  {
    Action->SetParticle(newValue);
  }
  else if (command == EnergyCmd)
  {
    Action->SetEnergy(EnergyCmd->GetNewDoubleValue(newValue));
  }
  else if (command == BeamSeedCmd)
  {
    Action->SetBeamSeed(BeamSeedCmd->GetNewIntValue(newValue));
  }
  else if (command == BeamProfileFileCmd)
  {
    Action->SetBeamProfileFile(newValue);
  }
  else if (command == TiltedBeamCmd) {
	Action->SetTiltAngle(TiltedBeamCmd->GetNewDoubleValue(newValue));
  }
  else if (command == RotCenterCmd) {
	Action->SetRotCenter(RotCenterCmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
