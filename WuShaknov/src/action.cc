//source for the initialization of user def routines

#include "action.hh"

MyActInit::MyActInit() {}
MyActInit::~MyActInit() {}

void MyActInit::Build() const {
	
	MyPrimGen* gen = new MyPrimGen();
	SetUserAction(gen);
	MyRunAction* runAction = new MyRunAction();
	SetUserAction(runAction);
	MyEventAction* eventAction = new MyEventAction(runAction);
	SetUserAction(eventAction);
	MySteppingAction* steppingAction = new MySteppingAction(eventAction);
	SetUserAction(steppingAction);

	G4cout << "RunAction set successfully!" << G4endl;
}

void MyActInit::BuildForMaster() const {

	MyRunAction* runAction = new MyRunAction();
	SetUserAction(runAction);
}