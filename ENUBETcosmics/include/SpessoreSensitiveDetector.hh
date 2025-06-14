#ifndef SpessoreSensitiveDetector_h
#define SpessoreSensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"

class G4Step;
class G4TouchableHistory;
class TFile;
class TH1D;



class SpessoreSensitiveDetector : public G4VSensitiveDetector
{
public:

  SpessoreSensitiveDetector(G4String SDname);
  virtual ~SpessoreSensitiveDetector();
  virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);

private:
    TFile*    m_ROOT_file;
};

#endif
