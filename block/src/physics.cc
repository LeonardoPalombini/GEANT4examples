//header for physics packages loading: EM physics, low E specific

#include "physics.hh"

MyPhysList::MyPhysList() {
	RegisterPhysics(new G4EmLivermorePhysics());
}
MyPhysList::~MyPhysList() {}