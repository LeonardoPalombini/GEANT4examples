//header for the run action

#ifndef RUN_HH
#define RUN_HH

#include "G4UserRunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4RootAnalysisManager.hh"


class MyRunAction : public G4UserRunAction
{
public:
	MyRunAction();
	~MyRunAction();

	virtual void BeginOfRunAction(const G4Run*);			//user def action to repeat once in a run
	virtual void EndOfRunAction(const G4Run*);
};


#endif