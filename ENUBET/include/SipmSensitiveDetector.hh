#ifndef SipmSensitiveDetector_h
#define SipmSensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"

class G4Step;
class G4TouchableHistory;

class SipmSensitiveDetector : public G4VSensitiveDetector
{
public:

  SipmSensitiveDetector(G4String SDname);
  virtual ~SipmSensitiveDetector();
  virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
};

#endif

