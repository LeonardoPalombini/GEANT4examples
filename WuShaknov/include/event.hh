//header for single-event routine

#ifndef EVENT_HH
#define EVENT_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "G4AnalysisManager.hh"
#include "G4RootAnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include "run.hh"

class MyEventAction :public G4UserEventAction {
public:
	MyEventAction(MyRunAction*);
	~MyEventAction();

	virtual void BeginOfEventAction(const G4Event*);
	virtual void EndOfEventAction(const G4Event*);

	void AddEdepR(G4double edep) { fEdepR += edep; }
	void AddEdepF(G4double edep) { fEdepF += edep; }


	void SetMonPosition(G4ThreeVector vec) { fmonPosition = vec; }

	void SetStepNo(G4int stpno) { fStepNo=stpno; }

	void SetID(G4int id) { fID = id; }

	void AddScattStepR() { fScatStepsR++; }
	void AddScattStepF() { fScatStepsF++; }

private:
	G4double fEdepR;
	G4double fEdepF;


	G4ThreeVector fmonPosition;

	G4int fStepNo;
	G4int fID;

	G4int fScatStepsR;
	G4int fScatStepsF;

};

#endif