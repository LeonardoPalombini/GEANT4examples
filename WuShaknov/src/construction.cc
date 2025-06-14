//source for det construction

#include "construction.hh"

MyDetConstr::MyDetConstr(double theta)
{
    thetadet = theta;
}

MyDetConstr::~MyDetConstr()
{}

G4VPhysicalVolume* MyDetConstr::Construct()
{
    G4double aluD = 3.6;    //alu cyl diameter cm
    G4double aluZ = 1.2;    //alu cyl height cm
    G4double aluPosX = 0.;  //alu cyl pos X cm
    G4double aluPosY = 0.;  //alu cyl pos Y cm
    G4double aluPosZF = 10.+1.9;    //fixed side alu cyl pos Z cm
    G4double aluPosZR = -10.-1.9;   //rotating side alu cyl pos Z cm

    G4double scintD = 5.1;   //scint diameter cm
    G4double scintZ = 5.1;  //scint height cm
    G4double scintRad = 10.5+2.55; //scint distance from axis cm
    G4double scintThR = thetadet; //rotating scint angle deg
    G4double scintPosZF = 10. + 2.55;  //fixed scint pos Z cm
    G4double scintPosZR = -10. - 2.55; //rotating scint pos Z cm

    G4double collX = 15.;   //collim X cm
    G4double collY = 21.;   //collim Y cm
    G4double collZ = 20.;   //collim Z cm
    G4double collL = 1.;    //collim hole D cm


    G4double scintPosXF = 0.;
    G4double scintPosYF = scintRad;
    G4double scintPosXR = scintRad * sin(scintThR * 3.1416 / 180.);
    G4double scintPosYR = scintRad * cos(scintThR * 3.1416 / 180.);
    


    G4NistManager* nist = G4NistManager::Instance();

    G4Material* worldMat = nist->FindOrBuildMaterial("G4_AIR");

    G4Box* solidWorld = new G4Box("solidWorld", .25 * m, .25 * m, .25 * m); //half-lengths
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
    G4VPhysicalVolume* physWorld = new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true); //param: rotation, position(x,y,z), obj, name, parent, bool op, copy ident, overlap check

    //alu block
    G4Material* alu = nist->FindOrBuildMaterial("G4_Al");

    G4RotationMatrix* rotAluR = new G4RotationMatrix();
    rotAluR->rotateX(90 * deg);
    rotAluR->rotateY(-scintThR * deg);
    G4RotationMatrix* rotAluF = new G4RotationMatrix();
    rotAluF->rotateX(90 * deg);

    G4Tubs* aluBlock = new G4Tubs("aluBlock", 0. * cm, 0.5 * aluD * cm, 0.5 * aluZ * cm, 0., 360. * deg); //inner R, outer R, half Z, start theta, end theta
    G4LogicalVolume* logicBlock = new G4LogicalVolume(aluBlock, alu, "logicBlock");
    G4VPhysicalVolume* physBlockF = new G4PVPlacement(rotAluF, G4ThreeVector(aluPosX, aluPosY, aluPosZF * cm), logicBlock, "physBlockF", logicWorld, false, 0, true);
    G4VPhysicalVolume* physBlockR = new G4PVPlacement(rotAluR, G4ThreeVector(aluPosX, aluPosY, aluPosZR * cm), logicBlock, "physBlockR", logicWorld, false, 0, true);

    //scintillator
    G4Element* lanthanum = nist->FindOrBuildElement("La");
    G4Element* bromine = nist->FindOrBuildElement("Br");
    G4Element* cerium = nist->FindOrBuildElement("Ce");

    G4Material* LaBr3 = new G4Material("LaBr3",5.1*g/cm3,3);
    LaBr3->AddElement(lanthanum, 19);//19
    LaBr3->AddElement(bromine, 57);//57
    LaBr3->AddElement(cerium, 1);   //cerium-activated, doping 0.5% in mass

    //collimator
    G4Material* lead = nist->FindOrBuildMaterial("G4_Pb");

    G4RotationMatrix* rotR= new G4RotationMatrix();
    rotR->rotateX(90 * deg);
    rotR->rotateY(-scintThR * deg);
    G4RotationMatrix* rotF = new G4RotationMatrix();
    rotF->rotateX(90 * deg);
    G4Tubs* scint = new G4Tubs("scint", 0. * cm, 0.5 * scintD * cm, 0.5 * scintZ * cm, 0., 360. * deg); //inner R, outer R, half Z, start theta, end theta
    logicScint = new G4LogicalVolume(scint, LaBr3, "logicScint");
    G4VPhysicalVolume* physScintF = new G4PVPlacement(rotF, G4ThreeVector(scintPosXF * cm, scintPosYF * cm, scintPosZF * cm), logicScint, "physScintF", logicWorld, false, 0, true);
    G4VPhysicalVolume* physScintR = new G4PVPlacement(rotR, G4ThreeVector(scintPosXR * cm, scintPosYR * cm, scintPosZR * cm), logicScint, "physScintR", logicWorld, false, 0, true);

    /*G4Tubs* wind = new G4Tubs("wind", 0. * cm, 0.5 * scintD * cm, 0.25 * mm, 0., 360. * deg); //inner R, outer R, half Z, start theta, end theta
    G4LogicalVolume* logicWind = new G4LogicalVolume(wind, alu, "logicWind");
    G4VPhysicalVolume* physWindF = new G4PVPlacement(rotF, G4ThreeVector(scintPosXF * cm, scintPosYF * cm - (0.25 * mm + 0.5 * scintZ*  cm), scintPosZF * cm), logicWind, "physWindF", logicWorld, false, 0, true);
    G4VPhysicalVolume* physWindR = new G4PVPlacement(rotR, G4ThreeVector(scintPosXR * cm, scintPosYR * cm - (0.25 * mm + 0.5 * scintZ * cm), scintPosZR * cm), logicWind, "physWindR", logicWorld, false, 0, true);
    */

    //source
    G4Material* plexi = nist->FindOrBuildMaterial("G4_PLEXIGLASS");

    G4Box* source = new G4Box("source",0.3*cm,0.3*cm,0.05*cm);
    G4LogicalVolume* logicSource = new G4LogicalVolume(source,plexi,"logicSource");
    G4VPhysicalVolume* physSource = new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), logicSource, "physSource", logicWorld, false, 0, true);

    
    
    G4Box* leadBlock = new G4Box("leadBlock", 0.5 * collX * cm, 0.5 * collY * cm, 0.5 * collZ * cm);
    G4Box* srcNeg = new G4Box("srcNeg", 0.35 * cm, 0.35 * cm, 0.1 * cm);
    G4Box* hole = new G4Box("hole", 0.5 * collL* cm, 0.5 * collL * cm, 0.5 * collZ * cm);
    G4SubtractionSolid* med = new G4SubtractionSolid("med",leadBlock,srcNeg);
    G4SubtractionSolid* collimator = new G4SubtractionSolid("collimator", med, hole);
    G4LogicalVolume* logicCollim = new G4LogicalVolume(collimator, lead, "logicCollim");
    G4VPhysicalVolume* physCollim = new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), logicCollim, "physCollim", logicWorld, false, 0, true);

    //monitor to study the collimator
    //G4Box* monitor = new G4Box("monitor",4*cm,4*cm,0.2*cm);
    //logicMon = new G4LogicalVolume(monitor,worldMat,"logicMon");
    //G4VPhysicalVolume* physMon = new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 10.25*cm),logicMon,"physMon",logicWorld, false, 0, true);
    
    physDets.id[0] = physScintR;
    physDets.id[1] = physScintF;

    getCollim = physCollim;
    //getMonitor = physMon;
    getSource = physSource;
    getScatt1 = physBlockF;
    getScatt2 = physBlockR;

    G4cout << "Setup built successfully!" << G4endl;
    return physWorld;
}


void MyDetConstr::ConstructSDandField() {

    MySensitiveDetector* sensScint = new MySensitiveDetector("sensScint");
    logicScint->SetSensitiveDetector(sensScint);

    //only for collimator study
    //MySensitiveDetector* sensMon = new MySensitiveDetector("sensMon");
    //logicMon->SetSensitiveDetector(sensMon);
    /*MySensitiveDetector* sensMoni = new MySensitiveDetector("sensMoni");
    logicScint->SetSensitiveDetector(sensMon);*/

}