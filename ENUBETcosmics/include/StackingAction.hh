#ifndef STACKINACTION_H
#define STACKINACTION_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class G4Track;

class StackingAction : public G4UserStackingAction {

public:
	StackingAction();
	virtual ~StackingAction();
	virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack );

};

#endif

