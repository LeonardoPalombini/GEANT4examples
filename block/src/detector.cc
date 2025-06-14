//source for sensitive detector action: hit energy and position

#include "detector.hh"
//#include <math.h>



MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name) {
}

MySensitiveDetector::~MySensitiveDetector() {}

G4bool MySensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory* R0hist) {


	return true;

}