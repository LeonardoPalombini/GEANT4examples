//header for sensitive detector action

#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4SystemOfUnits.hh"

#include "G4VSensitiveDetector.hh"

#include "G4AnalysisManager.hh"
#include "G4RootAnalysisManager.hh"

#include "G4RunManager.hh"

#include "G4PhysicsOrderedFreeVector.hh"




class MySensitiveDetector : public G4VSensitiveDetector
{
public:
	MySensitiveDetector(G4String);
	~MySensitiveDetector();

private:
	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

};




#endif
