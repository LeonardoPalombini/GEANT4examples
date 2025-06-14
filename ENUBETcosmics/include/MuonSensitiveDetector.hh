#ifndef MuonSensitiveDetector_h
#define MuonSensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"

class G4Step;
class G4TouchableHistory;
class TFile;
class TH1D;



class MuonSensitiveDetector : public G4VSensitiveDetector
{
public:

  MuonSensitiveDetector(G4String SDname);
  virtual ~MuonSensitiveDetector();
  virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);

private:
    TFile*    m_ROOT_file;
};

#endif
