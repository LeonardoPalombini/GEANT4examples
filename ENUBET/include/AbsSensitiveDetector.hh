#ifndef AbsSensitiveDetector_h
#define AbsSensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"

class G4Step;
class G4TouchableHistory;

class AbsSensitiveDetector : public G4VSensitiveDetector
{
public:

  AbsSensitiveDetector(G4String SDname);
  virtual ~AbsSensitiveDetector();
  virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
};

#endif

