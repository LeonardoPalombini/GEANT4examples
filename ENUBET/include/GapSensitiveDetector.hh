#ifndef GapSensitiveDetector_h
#define GapSensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"

class G4Step;
class G4TouchableHistory;
class G4HCofThisEvent;
class TFile;
class TH1D;



class GapSensitiveDetector : public G4VSensitiveDetector
{
public:

  GapSensitiveDetector(G4String SDname);
  ~GapSensitiveDetector();
  void Initialize(G4HCofThisEvent* hce);
  G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);

private:
    TFile*    m_ROOT_file;  

    std::map<G4int, std::vector<G4int>> gapCoord;
};

#endif
