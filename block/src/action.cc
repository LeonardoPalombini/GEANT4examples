//source for the initialization

#include "action.hh"

MyActInit::MyActInit(G4double en) {
	gammaEnergy = en;

}
MyActInit::~MyActInit() {}

void MyActInit::Build() const {
	MyPrimGen* gen = new MyPrimGen(gammaEnergy);
	SetUserAction(gen);

	MyRunAction* runAction = new MyRunAction();
	SetUserAction(runAction);

	MyEventAction* eventAction = new MyEventAction(runAction);
	SetUserAction(eventAction);

	MySteppingAction* steppingAction = new MySteppingAction(eventAction);
	SetUserAction(steppingAction);

}

void MyActInit::BuildForMaster() const {
	MyRunAction* runAction = new MyRunAction();
	SetUserAction(runAction);
}