//header for det construction

#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "detector.hh"

struct PhysDets {
    G4VPhysicalVolume* id;
};


class MyDetConstr : public G4VUserDetectorConstruction
{
public:
    MyDetConstr(double t);
    ~MyDetConstr();

    G4VPhysicalVolume* Construct();

    PhysDets GetPhysDetectors() const { return physDets; }

private:
    
    G4LogicalVolume* logicWall;
    PhysDets physDets;
    G4double Wthickness;

    virtual void ConstructSDandField();

};




#endif