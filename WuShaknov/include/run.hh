//header for the run action

#ifndef RUN_HH
#define RUN_HH

#include "G4UserRunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4RootAnalysisManager.hh"
#include "G4SystemOfUnits.hh"


class MyRunAction : public G4UserRunAction
{
public:
	MyRunAction();
	~MyRunAction();

	virtual void BeginOfRunAction(const G4Run*);			//user def action to repeat once in a run
	virtual void EndOfRunAction(const G4Run*);

	void addFEP() { iFEP = iFEP+1; };

	G4int iFEP;

private:

	
};


#endif