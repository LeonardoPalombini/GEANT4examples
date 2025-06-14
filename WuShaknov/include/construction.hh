//header for det construction

#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMAtrix.hh"
#include "G4SubtractionSolid.hh"

#include "detector.hh"

struct PhysDets {
    G4VPhysicalVolume* id[2];
};

class MyDetConstr : public G4VUserDetectorConstruction
{
public:
    MyDetConstr(double theta);
    ~MyDetConstr();

    G4VPhysicalVolume* Construct();
    PhysDets GetPhysDetectors() const { return physDets; }
    G4VPhysicalVolume* GetCollimator() const { return getCollim; }
    G4VPhysicalVolume* GetMonitor() const { return getMonitor; }
    G4VPhysicalVolume* GetSource() const { return getSource; }
    G4VPhysicalVolume* GetScatt1() const { return getScatt1; }
    G4VPhysicalVolume* GetScatt2() const { return getScatt2; }

private:
    G4double thetadet;
    G4LogicalVolume* logicScint;
    G4LogicalVolume* logicMon;
    virtual void ConstructSDandField();
    PhysDets physDets;
    G4VPhysicalVolume* getCollim;
    G4VPhysicalVolume* getMonitor;
    G4VPhysicalVolume* getSource;
    G4VPhysicalVolume* getScatt1;
    G4VPhysicalVolume* getScatt2;

};


#endif