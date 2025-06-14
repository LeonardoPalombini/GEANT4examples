//header for the initialization
#ifndef ACTION_HH
#define ACTION_HH

#include "G4VUserActionInitialization.hh"

#include "generator.hh"
#include "run.hh"
#include "event.hh"
#include "stepping.hh"

class MyActInit : public G4VUserActionInitialization
{
public:
	MyActInit(G4double en);
	~MyActInit();

	virtual void Build() const;
	virtual void BuildForMaster() const;

private:
	G4double gammaEnergy;

};

#endif