//source for det construction: alu cylinder, endcap detector for no-scattering events, 90° ring detector for transverse scattering events
#include "construction.hh"

MyDetConstr::MyDetConstr(double t)
{
    Wthickness = t;
}

MyDetConstr::~MyDetConstr()
{}

G4VPhysicalVolume* MyDetConstr::Construct()
{
    G4double scintD = 2.54 * 2.0;
    G4double scintZ = scintD;
    G4double scintPosX = 0.;
    G4double scintPosY = 0.;
    G4double scintPosZ = 120;
    G4double wallDist = 30.;
    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateX(0. * deg);

    
    G4NistManager* nist = G4NistManager::Instance();

    G4Material* worldMat = nist->FindOrBuildMaterial("G4_AIR");
    G4Box* solidWorld = new G4Box("solidWorld", 1. * m, 1.5 * m, 1.5 * m); //half-lengths
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
    G4VPhysicalVolume* physWorld = new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true); //param: rotation, position(x,y,z), obj, name, parent, bool op, copy ident, overlap check

    /*scintillator
    G4Element* lanthanum = nist->FindOrBuildElement("La");
    G4Element* bromine = nist->FindOrBuildElement("Br");
    G4Element* cerium = nist->FindOrBuildElement("Ce");

    G4Material* LaBr3 = new G4Material("LaBr3", 5.1 * g / cm3, 3);
    LaBr3->AddElement(lanthanum, 19);//19
    LaBr3->AddElement(bromine, 57);//57
    LaBr3->AddElement(cerium, 1);   //cerium-activated, doping 0.5% in masshysicalVolume* physRoundDet = new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 20. * cm), logicRoundDet, "physRoundDet", logicWorld, false, 0, true);
    
    G4Tubs* scint = new G4Tubs("scint", 0. * cm, 0.5 * scintD * cm, 0.5 * scintZ * cm, 0., 360. * deg); //inner R, outer R, half Z, start theta, end theta
    logicScint = new G4LogicalVolume(scint, LaBr3, "logicScint");
    G4VPhysicalVolume* physScint = new G4PVPlacement(rot, G4ThreeVector(scintPosX * cm, scintPosY * cm, scintPosZ * cm), logicScint, "physScintF", logicWorld, false, 0, true);

    physDets.id = physScint;
    

    //casing
    G4Material* Al = nist->FindOrBuildMaterial("G4_W");
    
    G4Tubs* casExt = new G4Tubs("casExt", 0. * cm, 0.5 * scintD * cm + 0.05*cm, 0.5 * scintZ * cm + 0.05*cm, 0., 360. * deg);
    G4SubtractionSolid* casing = new G4SubtractionSolid("casing", casExt, scint);
    G4LogicalVolume* logicCasing = new G4LogicalVolume(casing, Al, "logicCasing");
    G4VPhysicalVolume* physCasing = new G4PVPlacement(nullptr, G4ThreeVector(scintPosX * cm, scintPosY * cm, scintPosZ * cm), logicCasing, "physCasing", logicWorld, false, 0, true);
    */
    //chamber wall
    G4Material* W = nist->FindOrBuildMaterial("G4_W");
    G4Box* wall = new G4Box("wall", 15*cm, 15*cm, Wthickness*mm);
    logicWall = new G4LogicalVolume(wall, W, "logicWall");
    G4VPhysicalVolume* physWall = new G4PVPlacement(nullptr, G4ThreeVector(scintPosX * cm, scintPosY * cm, (scintPosZ-wallDist) * cm), logicWall, "physWall", logicWorld, false, 0, true);
    physDets.id = physWall;

    return physWorld;
}


void MyDetConstr::ConstructSDandField() {

    MySensitiveDetector* sensDet = new MySensitiveDetector("sensEndDet");
    logicWall->SetSensitiveDetector(sensDet);

}