
#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorMessenger : public G4UImessenger
{
public:
  PrimaryGeneratorMessenger(PrimaryGeneratorAction *);
  virtual ~PrimaryGeneratorMessenger();

  void SetNewValue(G4UIcommand *, G4String);

private:
  PrimaryGeneratorAction *Action;
  G4UIdirectory *gunDir;
  G4UIcmdWithAString *ProductionCmd;
  G4UIcmdWithAString *BeamTypeCmd;
  G4UIcmdWithAString *BeamProfileFileCmd;
  G4UIcmdWithAString *ParticleCmd;
  G4UIcmdWithAnInteger *BeamSeedCmd;
  G4UIcmdWithADouble *EnergyCmd;
	
  G4UIcmdWithADouble *TiltedBeamCmd;
  G4UIcmdWithADouble *RotCenterCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
