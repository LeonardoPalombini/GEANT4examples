//header for physics packages loading

#ifndef PHYSICS_HH
#define PHYSICS_HH

#include "G4VModularPhysicsList.hh"
#include "G4EmLivermorePhysics.hh"


class MyPhysList : public G4VModularPhysicsList
{
public:
	MyPhysList();
	~MyPhysList();
};



#endif