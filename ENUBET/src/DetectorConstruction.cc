#include "DetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"

#include "G4MaterialTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

#include "SipmSensitiveDetector.hh"
#include "AbsSensitiveDetector.hh"
#include "GapSensitiveDetector.hh"
#include "MuonSensitiveDetector.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include <iomanip>
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4Trd.hh"

#include "G4Version.hh"

#include "TRandom3.h"

#include "G4Types.hh"
#include "Randomize.hh"

#include "G4SDManager.hh"



G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0;

DetectorConstruction::DetectorConstruction()
:worldMaterial(0),solidWorld(0),logicWorld(0),physiWorld(0)
{

  nbOfBoxes = 1;
  
//the beam is along the X axis. Fibers go along Z axis 

  //Absorbers 
  absX = 15.*mm;
  absZ = 30*3*mm; 
  absY = 30*mm; 

  //Scintillators 
  gapX = 7*mm;
  gapY = 30 *mm; 
  gapZ = 30 *mm;

  //mother Volume for rotation 
  motherX = 1000*cm; 
  motherY = 1000*cm;
  motherZ = 1000*cm;
  
  

  //fibers
  //fiberDiameter = 1*mm; // ????????????????? or 1*mm ??
  fiberDiameter = 1*mm; // ????????????????? or 1*mm ??
  fiberDistance = 15*mm;
  fiberLenght3 = 45*cm+90*mm;
  fiberLenght2 = 45*cm+120*mm;
  fiberLenght1 = 45*cm+150*mm;
  fiberLenght0 = 45*cm + 180*mm;
  sensPlaneThickness = 5.0*mm;
  AirGapThickness = 10.0*mm;

  grooveDiameter = 1.1*mm; 
  
  grooveX = absX/2.+fiberDiameter/2.+0.5*mm;
  grooveZ = 30*cm;  
  grooveY = fiberDistance+6.4*mm; 
  //polietilene borato
  poliX = absX+gapX; 
  poliY = 30*mm; 
  poliZ = 30*cm; 

  nbOfXCalo=1;
  nbOfYCalo=2;
  nbOfZCalo=10;

  gapNum = 1;
  absNum=1;
  poliNum=1;
  fibNum=10000;

  ComputeCalorParameters();
  DefineMaterials();

  fCheckOverlaps=false;
  //fCheckOverlaps=false;

  //******* IRON ABSORBER DIMENSIONS ********
  arc_angle = pi/2.*radian;
  start_phi = pi/4.*radian;
  D_Ra = 11*cm;
  zAbs = 1.4*cm;
  layers_dist = 0.85*cm;
  nLayers = 75;
  demo_height = 36.5*cm;
    
  //****** BORATED POLYETHYLENE DIMENSIONS *******
  D_Rb = 30*cm;
  zBorPoly = 2.25*cm;
    
  //****** LAYERS INIZIALIZATION *********
  R_layers = 1;
  Phi_layers = 10;
  
  // ************** Muon Catcher ***************
  
  //muoncatcher fittizio: dimensioni

  MuonBoxX1 = 46.0*cm;
  MuonBoxY1 = 20.9*cm;
  MuonBoxZ1 = 4.2*cm;

  MuonBoxX2 = 27.2*cm;
  MuonBoxY2 = 12.0*cm;
  MuonBoxZ2 = 3.0*cm;
  
  MuonActiveVolX1 = 19.3*cm;
  MuonActiveVolY1 = 12.5*cm;
  MuonActiveVolZ1 = 3.0*cm;

  MuonActiveVolX2 = 10.0*cm;
  MuonActiveVolY2 = 10.0*cm;
  MuonActiveVolZ2 = 1.3*cm;
  
  // reasonable values
  BrickSizeX = 37*cm;
  BrickSizeY = 21*cm;
  BrickSizeZ = 6*cm;

  /*
  //spessore tra i due muon catcher
  FeX = 200*mm;
  FeY = 300*mm;//123*mm;
  FeZ = 300*mm;//100*mm;
  */
  
  //CLHEP::RandGauss rnd = CLHEP::RandGauss( *(CLHEP::HepRandom::getTheEngine()) , 0.0 , 1.0 );
  //rnd = CLHEP::RandGauss( *(CLHEP::HepRandom::getTheEngine()) , 0.0 , 1.0 );
  
  // sigma for randomization
  sigma = 0.12*mm;
  THICK_OFFSET = 0*mm;
  //G4double OffsetForConstruction = 50*cm;
  
  
}


DetectorConstruction::~DetectorConstruction() {}

void DetectorConstruction::DefineMaterials() {

  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* vacuum = nistManager->FindOrBuildMaterial("G4_Galactic");
  G4Material* air = nistManager->FindOrBuildMaterial("G4_AIR");
  G4Material* iron = nistManager->FindOrBuildMaterial("G4_Fe");
  G4Material* aluminum = nistManager->FindOrBuildMaterial("G4_Al");
  G4Material* silicon = nistManager->FindOrBuildMaterial("G4_Si");
  G4Material* lead = nistManager->FindOrBuildMaterial("G4_Pb");
  G4Material* concrete = nistManager->FindOrBuildMaterial("G4_CONCRETE");

  G4Element* H  = nistManager->FindOrBuildElement(1);
  G4Element* C  = nistManager->FindOrBuildElement(6);
  G4Element* O  = nistManager->FindOrBuildElement(8);
  G4Element* Si = nistManager->FindOrBuildElement(14);
  G4Element* B = nistManager->FindOrBuildElement(5);
  /*
  G4Material* scintillator = new G4Material("Scintillator",1.05*g/cm3, 2);
  scintillator->AddElement(C,8);
  scintillator->AddElement(H,8);
  scintillator->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
  */

  // Borated PolyEthylene - ENUBET - 5% Boron                                                                                                                                                                                                                          
  //density=1.04*g/cm3;
  G4Material* fBPE = new G4Material("BPE", 1.04*g/cm3, 3);
  fBPE->AddElement(H,88);
  fBPE->AddElement(C, 44);
  fBPE->AddElement(B, 3);


  G4Material* scintillator = new G4Material("Scintillator",1.03*g/cm3, 4);
  scintillator->AddElement(C,7);
  scintillator->AddElement(H,8);
  scintillator->AddElement(O,1);
  scintillator->AddElement(Si,1);
  scintillator->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  G4Material* fiberCore = nistManager->FindOrBuildMaterial("G4_POLYSTYRENE");
  G4Material* fiberClad1 = nistManager->FindOrBuildMaterial("G4_PLEXIGLASS");
  G4Material* fiberClad2 = nistManager->FindOrBuildMaterial("G4_POLYVINYLIDENE_FLUORIDE");

  nistManager->SetVerbose(0);


  //Material properties tables

  const G4int numPhotonEnergy=30;

  G4double photonEnergy[numPhotonEnergy]={2.066*eV,2.101*eV,2.137*eV,2.175*eV,2.214*eV,
                                          2.254*eV,2.296*eV,2.339*eV,2.384*eV,2.431*eV,
                                          2.455*eV,2.479*eV,2.530*eV,2.556*eV,2.583*eV,
                                          2.610*eV,2.637*eV,2.695*eV,2.724*eV,2.755*eV,
                                          2.786*eV,2.818*eV,2.850*eV,2.883*eV,2.917*eV,
                                          2.952*eV,2.987*eV,3.024*eV,3.099*eV,3.179*eV};

  G4double vacuumRInd[numPhotonEnergy]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
                                        1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
                                        1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};

  G4MaterialPropertiesTable* vacuumTable = new G4MaterialPropertiesTable();
  vacuumTable->AddProperty("RINDEX", photonEnergy, vacuumRInd,numPhotonEnergy);
  vacuum->SetMaterialPropertiesTable(vacuumTable);
  air->SetMaterialPropertiesTable(vacuumTable);

  G4double scintiRInd[numPhotonEnergy]={1.58,1.58,1.58,1.58,1.58,1.58,
                                        1.58,1.58,1.58,1.58,1.58,1.58,
                                        1.58,1.58,1.58,1.58,1.58,1.58,
                                        1.58,1.58,1.58,1.58,1.58,1.58,
                                        1.58,1.58,1.58,1.58,1.58,1.58};

  G4double scintiAbs[numPhotonEnergy]={160.*cm,160.*cm,160.*cm,160.*cm,160.*cm,
                                       160.*cm,160.*cm,160.*cm,160.*cm,160.*cm,
                                       160.*cm,160.*cm,160.*cm,160.*cm,160.*cm,
                                       160.*cm,160.*cm,160.*cm,160.*cm,160.*cm,
                                       160.*cm,160.*cm,160.*cm,160.*cm,160.*cm,
                                       160.*cm,160.*cm,160.*cm,160.*cm,160.*cm};

  G4double scintiFast[numPhotonEnergy]={0.0,0.0,0.0,0.0,0.0,0.0,
                                        0.0,0.0,0.0,0.0,5.0,10.0,
                                        13.0,16.0,18.0,21.0,25.0,36.0,
                                        42.0,55.0,65.0,69.0,71.0,81.0,
                                        100.0,94.0,71.0,25.0,5.0,0.0};

  G4MaterialPropertiesTable* scintiTable = new G4MaterialPropertiesTable();
  scintiTable->AddProperty("RINDEX",photonEnergy,scintiRInd,numPhotonEnergy);
  scintiTable->AddProperty("ABSLENGTH",photonEnergy,scintiAbs,numPhotonEnergy);
  // TODO: A.B. optical physics managed differently in G4-v11
  #if(G4VERSION_NUMBER >= 1100)
    scintiTable->AddProperty("SCINTILLATIONCOMPONENT1",photonEnergy, scintiFast,numPhotonEnergy);
  #else
    scintiTable->AddProperty("FASTCOMPONENT",photonEnergy, scintiFast,numPhotonEnergy);
  #endif
  scintiTable->AddConstProperty("SCINTILLATIONYIELD",10/keV);
  scintiTable->AddConstProperty("RESOLUTIONSCALE",1.0);
  // TODO: A.B. optical physics managed differently in G4-v11
  #if(G4VERSION_NUMBER >= 1100)
    scintiTable->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.4*ns);
  #else
    scintiTable->AddConstProperty("FASTTIMECONSTANT", 2.4*ns);
  #endif
  scintillator->SetMaterialPropertiesTable(scintiTable);


  G4double fiberCoreRInd[numPhotonEnergy]={1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.6,
                                           1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.6,
                                           1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.6};

  G4double fiberCoreAbs[numPhotonEnergy]={3.0*m,3.0*m,3.0*m,3.0*m,3.0*m,3.0*m,
                                          3.0*m,3.0*m,3.0*m,3.0*m,3.0*m,3.0*m,
                                          175.0*cm,33.0*cm,21.0*cm,11.0*cm,5.0*cm,0.28*cm,
                                          0.27*cm,0.25*cm,0.22*cm,0.19*cm,0.13*cm,0.08*cm,
                                          0.06*mm,0.05*cm,0.03*cm,0.02*cm,0.03*cm,0.05*cm};

  G4double fiberCoreEmission[numPhotonEnergy]={5.0,9.0,10.0,14.0,19.0,26.0,
                                               35.0,49.0,62.0,82.0,90.0,97.0,
                                               100.0,94.0,70.0,45.0,12.0,0.0,
                                               0.0,0.0,0.0,0.0,0.0,0.0,
                                               0.0,0.0,0.0,0.0,0.0,0.0};

  G4MaterialPropertiesTable* fiberCoreTable = new G4MaterialPropertiesTable();
  fiberCoreTable->AddProperty("RINDEX",photonEnergy,fiberCoreRInd,numPhotonEnergy);
  fiberCoreTable->AddProperty("WLSABSLENGTH",photonEnergy,fiberCoreAbs,numPhotonEnergy);
  fiberCoreTable->AddProperty("WLSCOMPONENT",photonEnergy, fiberCoreEmission,numPhotonEnergy);
  fiberCoreTable->AddConstProperty("WLSTIMECONSTANT",2.7*ns);
  fiberCore->SetMaterialPropertiesTable(fiberCoreTable);

  G4double fiberCladAbs[numPhotonEnergy]={3.0*m,3.0*m,3.0*m,3.0*m,3.0*m,3.0*m,
                                          3.0*m,3.0*m,3.0*m,3.0*m,3.0*m,3.0*m,
                                          3.0*m,3.0*m,3.0*m,3.0*m,3.0*m,3.0*m,
                                          3.0*m,3.0*m,3.0*m,3.0*m,3.0*m,3.0*m,
                                          3.0*m,3.0*m,3.0*m,3.0*m,3.0*m,3.0*m};

  G4double fiberClad1RInd[numPhotonEnergy]={1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,
                                            1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,
                                            1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49};

  G4MaterialPropertiesTable* fiberClad1Table = new G4MaterialPropertiesTable();
  fiberClad1Table->AddProperty("RINDEX",photonEnergy,fiberClad1RInd,numPhotonEnergy);
  fiberClad1Table->AddProperty("ABSLENGTH",photonEnergy,fiberCladAbs,numPhotonEnergy);
  fiberClad1->SetMaterialPropertiesTable(fiberClad1Table);


  G4double fiberClad2RInd[numPhotonEnergy]={1.42,1.42,1.42,1.42,1.42,1.42,1.42,1.42,1.42,1.42,
                                            1.42,1.42,1.42,1.42,1.42,1.42,1.42,1.42,1.42,1.42,
                                            1.42,1.42,1.42,1.42,1.42,1.42,1.42,1.42,1.42,1.42};

  G4MaterialPropertiesTable* fiberClad2Table = new G4MaterialPropertiesTable();
  fiberClad2Table->AddProperty("RINDEX",photonEnergy,fiberClad2RInd,numPhotonEnergy);
  fiberClad2Table->AddProperty("ABSLENGTH",photonEnergy,fiberCladAbs,numPhotonEnergy);
  fiberClad2->SetMaterialPropertiesTable(fiberClad2Table);


  G4double siliconRInd[numPhotonEnergy]={1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,
                                         1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,
                                         1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4};

  G4MaterialPropertiesTable* siliconTable = new G4MaterialPropertiesTable();
  siliconTable->AddProperty("RINDEX",photonEnergy,siliconRInd,numPhotonEnergy);
  silicon->SetMaterialPropertiesTable(siliconTable);

  worldMaterial = air;//air
  absorMaterial = iron;//iron//lead
  gapMaterial = scintillator;
  boxMaterial = aluminum;
  boxMaterial1 = air;
  fiberMaterial = fiberCore;
  cladding1Material = fiberClad1;
  cladding2Material = fiberClad2;
  sipmMaterial = silicon;
  sensorPlaneMaterial = scintillator;//aluminum
  airCushionMaterial= air;
  brickMaterial = concrete;
  poliMaterial = fBPE;
    
  ironAbsMaterial = iron;
  borpolyMaterial = fBPE;
  
  MCMaterial = scintillator;//muoncatcher

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

void DetectorConstruction::ComputeCalorParameters() {

  calorX=0.;
  calorZ=0.;
  nbOfPlanes=0;

  calorX = nbOfAbsor*absX + nbOfGaps*gapX;
  calorZ = nbOfAbsor*absZ;
  totalCalorThickness = calorX+sensPlaneThickness+ AirGapThickness;
  nbOfPlanes = nbOfAbsor+nbOfGaps;

  worldSizeX = 10000*cm;
  worldSizeYZ = 10000*cm;

  fibersPerRow = G4int(sqrt(nbOfFibers));
  fiberStart = (absY-(fibersPerRow-1)*fiberDistance)/2.;

}

void DetectorConstruction::ComputeDemonstratorParameters() {
	
	// This function compute the demonstrator components
	
	R_t0 = 97*cm; // internal radius of t0-layer
	th_tile = 1.77*deg; // angle of a phi-sector of demonstrator
	D_Rg = 0.2*cm; // gap between absorbers and Borated polyethylene (BPE)
	N_phi_instr = 10; // number of instrumented phi-sectors
	N_z_instr = 8; // number of instrumented z layer (LCM units)
	N_LCM = 240; // total number of LCMs
	N_tile_LCM = 1200; // total number of tiles in the calorimeter
	N_tile_t0 = 160; // total number of tiles in t0-layer
	// tiles are trapezoidal
	
	h = 3*cm; // height
    	tile_thickness_nom = 6.7*mm;
    	
    	tile_thickness_eff = tile_thickness_nom + THICK_OFFSET;
	// In the following maj_base and min_base suffixes stand for major and minor base of the trapezoid
	// convert th_tile from deg to rad
	th_tile =  th_tile *rad;
	
	// t0
	t0_min_base = R_t0 * th_tile;
	t0_maj_base = ( R_t0 + h ) * th_tile;
	
	// t1 and t4
	t1_min_base = ( R_t0 + h ) * th_tile;
	t1_maj_base = ( R_t0 + 2. * h ) * th_tile;
	t4_min_base = t1_min_base;
	t4_maj_base = t1_maj_base;
	// t2 and t5
	t2_min_base = ( R_t0 +  2. * h ) * th_tile;
	t2_maj_base = ( R_t0 + 3. * h ) * th_tile;
	t5_min_base = t2_min_base;
	t5_maj_base = t2_maj_base;
	// t3 and t6
	t3_min_base = ( R_t0 +  3. * h ) * th_tile;
	t3_maj_base = ( R_t0 + 4. * h ) * th_tile;
	t6_min_base = t3_min_base;
	t6_maj_base = t3_maj_base;
	
	// absorbers	
	Rint_abs = R_t0 + h; // radius of inner arch
	arch_angle_abs = 90*deg;
	D_Ra = 11*cm;
	abs_thickness = 1.4*cm; // thickness of absorbers
	abs_distance = 0.85*cm; // distance among absorbers
	Ntot_abs = 75; // total number of absorbers
	demo_height = 36.5*cm; //height between ground and centre of the inner arc of iron absorber (along y-axis)
	
	// BPE
	Rint_BPE = R_t0 + h + D_Ra + D_Rg; // radius of inner arch
	arch_angle_BPE = 90*deg;
	D_Rb = 30*cm;
	BPE_thickness = 2.25*cm; // thickness of BPE
	Ntot_BPE = 75; // total number of BPE
}


G4VPhysicalVolume* DetectorConstruction::Construct() {

  solidWorld = new G4Box("solidWorld",
                         worldSizeX/2.,worldSizeYZ/2.,worldSizeYZ/2.);

  logicWorld = new G4LogicalVolume(solidWorld,
                                   worldMaterial,
                                   "logicWorld");
  //G4RotationMatrix* rmW = new G4RotationMatrix();
  //rmW->rotateY(90*deg);
  physiWorld = new G4PVPlacement(0,
  				                       G4ThreeVector(0.,0.,0.),
                                 logicWorld,
                                 "physiWorld",
                                 0,
                                 false,
                                 0);

  // TODO: A.B. managed differently in G4-v11
  #if(G4VERSION_NUMBER >= 1100)
    logicWorld->SetVisAttributes (G4VisAttributes::GetInvisible());
  #else
    logicWorld->SetVisAttributes (G4VisAttributes::Invisible);
  #endif

  ComputeDemonstratorParameters();
  //ConstructCalorimeter();
  ConstructTiles();
    
  ConstructIronAbs_BorPoly();
  ConstructMuCatcher();

  return physiWorld;
}


/*

void DetectorConstruction::ConstructCalorimeter() {

  

  G4VisAttributes* poliAttr= new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  G4VisAttributes* motherAttr= new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  G4VisAttributes* gapAttr= new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  G4VisAttributes* absAttr= new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  G4VisAttributes* fiberAttr= new G4VisAttributes(G4Colour(0.0,0.8,0.0));
  G4VisAttributes* clad1Attr= new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  G4VisAttributes* clad2Attr= new G4VisAttributes(G4Colour(0.0,0.5,0.0));
  G4VisAttributes* sipmAttr= new G4VisAttributes(G4Colour(0.8,0.8,0.8));
  G4VisAttributes* sPlaneAttr= new G4VisAttributes(G4Colour(0.1,0.6,0.1));

  clad1Attr->SetVisibility(true);
  clad1Attr->SetForceSolid(true);

  clad2Attr->SetVisibility(true);
  clad2Attr->SetForceSolid(true);

  absAttr->SetForceSolid(true);
  gapAttr->SetVisibility(true);
  poliAttr->SetVisibility(true);

  motherAttr->SetVisibility(false);

  fiberAttr->SetVisibility(true);
  fiberAttr->SetForceSolid(true);

  sipmAttr->SetForceSolid(true);

////////////SOLID VOLUMES//////////////////////////////
//mother volume
  //G4Box *motherBox = new G4Box("Mother",motherX/2.,motherY/2.,motherZ/2.);
//assorbitori
  //G4Box *absBox = new G4Box("Absorber",absX/2.,absY/2.,absZ/2.);

  
  
  
  // Quadruplet Builder
  
  //tiles are trapezoidal, grooves are rectangular
  
  
  // tile 0
  G4Trd* tile0 = new G4Trd("tile 0", tile_thickness_eff/2., tile_thickness_eff/2., t0_min_base/2., t0_maj_base/2., h/2.); 
  // G4Box *groove = new G4Box("innerfiber", grooveDiameter/2., grooveDiameter/2., h/2.);
  G4Box *groove = new G4Box("innerfiber", grooveDiameter/2. + THICK_OFFSET/4., grooveDiameter/2., h/2.); // asymmetrical groove !! THICKNESS OFFSET

  G4SubtractionSolid *hollowGroove_t0 = new G4SubtractionSolid("hollow groove of tile0", tile0, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  G4SubtractionSolid *t0 = new G4SubtractionSolid("t0", hollowGroove_t0, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  
  // tile 1
  G4Trd* tile1 = new G4Trd("tile 1", tile_thickness_eff/2., tile_thickness_eff/2., t1_min_base/2., t1_maj_base/2., h/2.);
  G4SubtractionSolid *hollowGroove_t1_1 = new G4SubtractionSolid("hollow groove of tile1 - 1", tile1, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  G4SubtractionSolid *hollowGroove_t1_2 = new G4SubtractionSolid("hollow groove of tile1 - 2", hollowGroove_t1_1, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  G4SubtractionSolid *hollowGroove_t1_3 = new G4SubtractionSolid("hollow groove of tile1 - 3", hollowGroove_t1_2, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2)/2., -fiberDistance/2.,0.));
  G4SubtractionSolid *t1 = new G4SubtractionSolid("t1", hollowGroove_t1_3, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  
  // tile 2
  G4Trd* tile2 = new G4Trd("tile 2", tile_thickness_eff/2., tile_thickness_eff/2., t2_min_base/2., t2_maj_base/2., h/2.);
  G4SubtractionSolid *hollowGroove_t2_1 = new G4SubtractionSolid("hollow groove of tile2 - 1", tile2, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  G4SubtractionSolid *hollowGroove_t2_2 = new G4SubtractionSolid("hollow groove of tile2 - 2", hollowGroove_t2_1, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  G4SubtractionSolid *hollowGroove_t2_3 = new G4SubtractionSolid("hollow groove of tile2 - 3", hollowGroove_t2_2, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  G4SubtractionSolid *hollowGroove_t2_4 = new G4SubtractionSolid("hollow groove of tile2 - 4", hollowGroove_t2_3, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  G4SubtractionSolid * hollowGroove_t2_5 = new G4SubtractionSolid("hollow groove of tile2 - 5", hollowGroove_t2_4, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. -2*mm,0.));
  G4SubtractionSolid *t2 = new G4SubtractionSolid("t2", hollowGroove_t2_5, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance -2*mm,0.));
  
  // tile 3
  G4Trd* tile3 = new G4Trd("tile 3", tile_thickness_eff/2., tile_thickness_eff/2., t3_min_base/2., t3_maj_base/2., h/2.);
  G4SubtractionSolid *hollowGroove_t3_1 = new G4SubtractionSolid("hollow groove of tile3 - 1", tile3, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  G4SubtractionSolid *hollowGroove_t3_2 = new G4SubtractionSolid("hollow groove of tile3 - 2", hollowGroove_t3_1, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  G4SubtractionSolid *hollowGroove_t3_3 = new G4SubtractionSolid("hollow groove of tile3 - 3", hollowGroove_t3_2, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  G4SubtractionSolid *hollowGroove_t3_4 = new G4SubtractionSolid("hollow groove of tile3 - 4", hollowGroove_t3_3, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  G4SubtractionSolid * hollowGroove_t3_5 = new G4SubtractionSolid("hollow groove of tile3 - 5", hollowGroove_t3_4, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. -2*mm,0.));
  G4SubtractionSolid * hollowGroove_t3_6 = new G4SubtractionSolid("hollow groove of tile3 - 6", hollowGroove_t3_5, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance -2*mm,0.));
  G4SubtractionSolid * hollowGroove_t3_7 = new G4SubtractionSolid("hollow groove of tile3 - 7", hollowGroove_t3_6, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + 2*mm,0.));
  G4SubtractionSolid *t3 = new G4SubtractionSolid("t3", hollowGroove_t3_7, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance + 2*mm,0.));
  
  
  
  
  // Triplet Builder
  
  //tiles are trapezoidal, grooves are rectangular
  
  
   
  
  // tile 4
  G4Trd* tile4 = new G4Trd("tile 4", tile_thickness_eff/2., tile_thickness_eff/2., t4_min_base/2., t4_maj_base/2., h/2.); 
  G4SubtractionSolid *hollowGroove_t4 = new G4SubtractionSolid("hollow groove of tile4", tile4, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  G4SubtractionSolid *t4 = new G4SubtractionSolid("t4", hollowGroove_t4, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  
  // tile 5
  G4Trd* tile5 = new G4Trd("tile 5", tile_thickness_eff/2., tile_thickness_eff/2., t5_min_base/2., t5_maj_base/2., h/2.); 
  G4SubtractionSolid *hollowGroove_t5_1 = new G4SubtractionSolid("hollow groove of tile4 - 1", tile5, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  G4SubtractionSolid *hollowGroove_t5_2 = new G4SubtractionSolid("hollow groove of tile4 -2", hollowGroove_t5_1, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  G4SubtractionSolid *hollowGroove_t5_3 = new G4SubtractionSolid("hollow groove of tile4 -3", hollowGroove_t5_2, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  G4SubtractionSolid *t5 = new G4SubtractionSolid("t5", hollowGroove_t5_3, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  
  // tile 6
  G4Trd* tile6 = new G4Trd("tile 6", tile_thickness_eff/2., tile_thickness_eff/2., t6_min_base/2., t6_maj_base/2., h/2.); 
  G4SubtractionSolid *hollowGroove_t6_1 = new G4SubtractionSolid("hollow groove of tile6 - 1", tile6, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  G4SubtractionSolid *hollowGroove_t6_2 = new G4SubtractionSolid("hollow groove of tile6 -2", hollowGroove_t6_1, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  G4SubtractionSolid *hollowGroove_t6_3 = new G4SubtractionSolid("hollow groove of tile6 -3", hollowGroove_t6_2, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  G4SubtractionSolid *hollowGroove_t6_4 = new G4SubtractionSolid("hollow groove of tile6 -4", hollowGroove_t6_3, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  G4SubtractionSolid *hollowGroove_t6_5 = new G4SubtractionSolid("hollow groove of tile6 -5", hollowGroove_t6_4, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance + 2*mm,0.));
  
  G4SubtractionSolid *t6 = new G4SubtractionSolid("t6", hollowGroove_t6_5, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + 2*mm,0.));
  
  */
  
  /*
//polietilene borato
  G4Box *outerPoli = new G4Box("tile 0",poliX/2., poliY/2., poliZ/2.);
  G4Box* groove_Poli = new G4Box("outer groove",grooveX/2., grooveY/2., grooveZ/2.); 
  //G4Tubs *innerfiber2 = new G4Tubs ("innerfiber2", 0, fiberDiameter/2.+0.5*mm, poliZ/2., 0, 2*M_PI*rad);
  G4SubtractionSolid *hollowPoli = new G4SubtractionSolid("Hollow Poli", outerPoli, groove_Poli, 0, G4ThreeVector(-(poliX/2.-0.55*mm),0.,0.));
  G4SubtractionSolid *poliFinal = new G4SubtractionSolid("Hollow Poli", hollowPoli, groove_Poli, 0, G4ThreeVector(+(poliX/2.-0.55*mm),0.,0.));
  
  
  // FIX ME: IT was like this before
  //polietilene borato
  G4Box *outerPoli = new G4Box("tile 0",poliX/2., poliY/2., poliZ/2.);
  G4Box* grooveBPE = new G4Box("outer groove",grooveX/2., grooveY/2., grooveZ/2.); 
  //G4Tubs *innerfiber2 = new G4Tubs ("innerfiber2", 0, fiberDiameter/2.+0.5*mm, poliZ/2., 0, 2*M_PI*rad);
   G4SubtractionSolid *hollowPoli = new G4SubtractionSolid("Hollow Poli", outerPoli, grooveBPE, 0, G4ThreeVector(-(poliX/2.-0.55*mm),0.,0.));
  G4SubtractionSolid *poliFinal = new G4SubtractionSolid("Hollow Poli", hollowPoli, grooveBPE, 0, G4ThreeVector(+(poliX/2.-0.55*mm),0.,0.));
  */
  
/*
  solidAbs = new G4Box("solidAbsorber",absX/2.,absY/2.,absZ/2.);
  solidGap = new G4Box("solidGap",gapX/2.,gapY/2.,gapZ/2.);
  solidPoli = new G4Box("solidPoli", poliX/2., poliY/2., poliZ/2.); 
  //scatola esterna
  solidMother = new G4Box("solidMother",motherX/2.,motherY/2.,motherZ/2.);
  

  //solidMother_quad = new G4Trd("solidMother_quad", tile_thickness_eff/2., tile_thickness_eff/2., t0_min_base/2., ( R_t0 + L ) * th_tile / 2. , L/2. );
  
  //solidMother_quad = new G4Trd("solidMother_quad", tile_thickness_eff/2., tile_thickness_eff/2., t0_min_base/2., t3_maj_base/2., 2*h );
  

  G4double sipmDimX=0.5*mm;
  G4double sipmDimY=fiberDiameter/2.;
  G4double sipmDimZ=fiberDiameter/2.;
  G4double airCushionDimX = sipmDimX/10.;

  solidSensorPlane = new G4Box("solidSensorPlane",sensPlaneThickness/2.,absY/2.,absZ/2.);
  solidSipm = new G4Box("solidSipm",sipmDimX,sipmDimY,sipmDimZ);
  solidAirCushion = new G4Box("solidAirCushion",airCushionDimX,sipmDimY,sipmDimZ);

  ////////////LOGIC VOLUMES//////////////////////////////
  logicMother = new G4LogicalVolume(solidMother,worldMaterial,"logicMother");
  //logicMother_quad = new G4LogicalVolume(solidMother_quad,worldMaterial,"logicMother_quad");
  logicAbs = new G4LogicalVolume(solidAbs,absorMaterial,"logicAbsorber");
  logict0 = new G4LogicalVolume(t0,gapMaterial,"logict0");
  logict1 = new G4LogicalVolume(t1,gapMaterial,"logict1");
  logict2 = new G4LogicalVolume(t2,gapMaterial,"logict2");
  logict3 = new G4LogicalVolume(t3,gapMaterial,"logict3");
  logict4 = new G4LogicalVolume(t4,gapMaterial,"logict4");
  logict5 = new G4LogicalVolume(t5,gapMaterial,"logict5");
  logict6 = new G4LogicalVolume(t6,gapMaterial,"logict6");
 // logicPoli = new G4LogicalVolume(poliFinal, poliMaterial, "poliFinal");

  //logicWLSAbs = ConstructWLSFiber(calorSizeYZ);
  // quadruplet
  logicWLSGap0 = ConstructWLSFiber(fiberLenght0);
  logicWLSGap1 = ConstructWLSFiber(fiberLenght1);
  logicWLSGap2 = ConstructWLSFiber(fiberLenght2);
  logicWLSGap3 = ConstructWLSFiber(fiberLenght3);
  // triplet
  logicWLSGap1_trip = ConstructWLSFiber(fiberLenght1);
  logicWLSGap2_trip = ConstructWLSFiber(fiberLenght2);
  logicWLSGap3_trip = ConstructWLSFiber(fiberLenght3);

  logicSensorPlane = new G4LogicalVolume(solidSensorPlane,sensorPlaneMaterial,"logicSensorPlane");
  logicSipm = new G4LogicalVolume(solidSipm,sipmMaterial,"logicSipm");
  logicAirCushion = new G4LogicalVolume(solidAirCushion,airCushionMaterial,"logicAirCushion");

  logicAbs->SetVisAttributes(absAttr);
  logict0->SetVisAttributes(gapAttr);
  logict1->SetVisAttributes(gapAttr);
  logict2->SetVisAttributes(gapAttr);
  logict3->SetVisAttributes(gapAttr);
  logict4->SetVisAttributes(gapAttr);
  logict5->SetVisAttributes(gapAttr);
  logict6->SetVisAttributes(gapAttr);
  //logicPoli->SetVisAttributes(poliAttr);
  logicMother->SetVisAttributes(motherAttr);

  logicWLSGap0->SetVisAttributes(fiberAttr);
  logicWLSGap1->SetVisAttributes(fiberAttr); 
  logicWLSGap2->SetVisAttributes(fiberAttr); 
  logicWLSGap3->SetVisAttributes(fiberAttr);
  // triplet
  logicWLSGap1_trip->SetVisAttributes(fiberAttr);
  logicWLSGap2_trip->SetVisAttributes(fiberAttr);
  logicWLSGap3_trip->SetVisAttributes(fiberAttr);
  logicFiber->SetVisAttributes(fiberAttr);
  logicClad1->SetVisAttributes(clad1Attr);
  logicClad2->SetVisAttributes(clad2Attr);
  logicSipm->SetVisAttributes(sipmAttr);
  logicSensorPlane->SetVisAttributes(sPlaneAttr);
  logicAirCushion->SetVisAttributes(clad1Attr);


  //////////////////PHYSICAL PLACEMENTS//////////////////////////////

  physiAirCushion =
      new G4PVPlacement(0,
                        G4ThreeVector(-sipmDimX+airCushionDimX,0,0),
                        logicAirCushion,
                        "physiAirCushion",
                        logicSipm,
                        false,
                        sipmNum+500,
                        fCheckOverlaps);

  G4double sipmXcenter = -sensPlaneThickness/2. + sipmDimX;

  for(G4int i=0;i<nbOfFibers;i++) {
    G4double ycenter=absY/2.-fiberStart-(i%fibersPerRow)*fiberDistance;
    G4double zcenter=absZ/2.-fiberStart-G4int(i/fibersPerRow)*fiberDistance;

    physiSipm =
        new G4PVPlacement(0,
                          G4ThreeVector(sipmXcenter,ycenter,zcenter),
                          logicSipm,
                          "physiSipm",
                          logicSensorPlane,
                          false,
                          sipmNum+i,
                          fCheckOverlaps);
  }


//G4RotationMatrix* rm = new G4RotationMatrix();
//rm->rotateX(90.*deg);
for(G4int nBoxloop=0;nBoxloop<nbOfBoxes;nBoxloop++)
{

// Rotation of the Mother volume

G4RotationMatrix* rmM = new G4RotationMatrix();

*/

/*rmM->rotateZ(90.*deg);
rmM->rotateX(-90.*deg);
rmM->rotateY(90.*deg);*/

/*
rmM->rotateY(-90.*deg); // just for quadruplet and triplet screenshootting



//rmM->rotateX(90.*deg);

//rmM->rotateZ(-90.*deg);
//rmM->rotateX(-90.*deg);
physiMother=
        new G4PVPlacement(rmM,
                          G4ThreeVector(0., 0., 0.),
                          logicMother,//logicWLSGap1
                          "physiMother",
                          logicWorld,
                          false,
                          0,
                      fCheckOverlaps);
 
 */ 
 /* 
G4RotationMatrix* rm = new G4RotationMatrix();
//rm->rotateY();
for(G4int xCalo=0;xCalo<nbOfXCalo;xCalo++) { //in this case we only have one but we leave it for extending to more in the future
for(G4int zCalo=0;zCalo<nbOfZCalo;zCalo++) {
for(G4int yCalo=0;yCalo<nbOfYCalo;yCalo++) {
  //nbOfFibers = 1;//5
  nbOfFibers = 1;//5
     for(G4int i=0;i<nbOfFibers;i++){
     
     
  //fibers coordinates and placement
    G4double fxcenter = -calorX/2.+gapX/2.+ (gapX+absX)*i + nBoxloop*totalCalorThickness;
    G4double fycenter = gapY/2. +yCalo*gapY-5.5*mm;
    G4double fzcenter = gapZ/2. +zCalo*gapZ;
    G4int Findex =zCalo*nbOfYCalo*nbOfFibers + yCalo*nbOfFibers + nBoxloop*nbOfZCalo*nbOfYCalo*nbOfFibers + i;
    
    
    */
   //G4RotationMatrix* rm = new G4RotationMatrix(0*deg, 90*deg, 90*deg);
	/*
    if(i<2){

    if (yCalo < 1){

      physiWLSGap[i]=
        new G4PVPlacement(0,
                          G4ThreeVector(fxcenter+gapX/2.-fiberDiameter/2.-0.2*mm,fycenter-2*mm,fzcenter+fiberLenght1/2.-gapZ/2.+0.1*mm -30*mm),
                          logicWLSGap1,//logicWLSGap1
                          "WLSfiber",
                          logicMother,
                          false,
                          fibNum+Findex,
                      fCheckOverlaps);

    physiWLSGap[i]=
        new G4PVPlacement(0,
                          G4ThreeVector(fxcenter+gapX/2.-fiberDiameter/2.-0.2*mm,fycenter -2*mm + fiberDistance,fzcenter+fiberLenght1/2.-gapZ/2.+0.1*mm -30*mm),
                          logicWLSGap1,//logicWLSGap1
                          "WLSfiber",
                          logicMother,
                          false,
                          fibNum+Findex,
                      fCheckOverlaps);

    }}*/
    
    
  /*  if(zCalo == 0){ 
  
    
    // boooo
    
    
    
    
    physiWLSGap[i]=
        new G4PVPlacement(0,
                          G4ThreeVector(fxcenter+gapX/2.-fiberDiameter/2.-0.2*mm,fycenter-2*mm,fzcenter+fiberLenght1/2.-gapZ/2.+0.1*mm -30*mm),
                          logicWLSGap1,//logicWLSGap1
                          "WLSfiber_quad",
                          logicMother,
                          false,
                          fibNum+Findex,
                      fCheckOverlaps);

    physiWLSGap[i]=
        new G4PVPlacement(0,
                          G4ThreeVector(fxcenter+gapX/2.-fiberDiameter/2.-0.2*mm,fycenter -2*mm + fiberDistance,fzcenter+fiberLenght1/2.-gapZ/2.+0.1*mm -30*mm),
                          logicWLSGap1,//logicWLSGap1
                          "WLSfiber_quad",
                          logicMother,
                          false,
                          fibNum+Findex,
                      fCheckOverlaps);
    	
    	// quadruplet
    	
        physiWLSGap[i]=
        new G4PVPlacement(0,
                          G4ThreeVector(fxcenter+fiberDiameter/2.-gapX/2.+0.2*mm,fycenter -2*mm,fzcenter+fiberLenght2/2.-gapZ/2.+0.1*mm),
                          logicWLSGap2,//logicWLSGap1
                          "WLSfiber_quad",
                          logicMother,
                          false,
                          fibNum+Findex,
                      fCheckOverlaps);

        physiWLSGap[i]=
        new G4PVPlacement(0,
                          G4ThreeVector(fxcenter+fiberDiameter/2.-gapX/2.+0.2*mm,fycenter -2*mm + fiberDistance,fzcenter+fiberLenght2/2. -gapZ/2.+ 0.1*mm),
                          logicWLSGap2,//logicWLSGap1
                          "WLSfiber_quad",
                          logicMother,
                          false,
                          fibNum+Findex,
                      fCheckOverlaps);
                      
       // tiplet             
       
       physiWLSGap_trip[i]=
        new G4PVPlacement(0,
                          G4ThreeVector(fxcenter+fiberDiameter/2.-gapX/2.+0.2*mm -OffsetForConstruction,fycenter -2*mm,fzcenter+fiberLenght2/2.-gapZ/2.+0.1*mm),
                          logicWLSGap1_trip,//logicWLSGap1
                          "WLSfiber_trip",
                          logicMother,
                          false,
                          fibNum+Findex,
                      fCheckOverlaps);

        physiWLSGap_trip[i]=
        new G4PVPlacement(0,
                          G4ThreeVector(fxcenter+fiberDiameter/2.-gapX/2.+0.2*mm -OffsetForConstruction,fycenter -2*mm + fiberDistance,fzcenter+fiberLenght2/2. -gapZ/2.+ 0.1*mm),
                          logicWLSGap1_trip,//logicWLSGap1
                          "WLSfiber_trip",
                          logicMother,
                          false,
                          fibNum+Findex,
                      fCheckOverlaps);
                      
                      
                      
                      
                 
      }

      if(zCalo == 1){
      
      	// quadruplet
        physiWLSGap[i]=
        new G4PVPlacement(0,
                          G4ThreeVector(fxcenter-fiberDiameter/2.+gapX/2.-0.2*mm,fycenter - 4*mm,fzcenter+fiberLenght3/2. -gapZ/2.+0.1*mm),
                          logicWLSGap3,//logicWLSGap1
                          "WLSfiber_quad",
                          logicMother,
                          false,
                          fibNum+Findex,
                      fCheckOverlaps);

        physiWLSGap[i]=
        new G4PVPlacement(0,
                          G4ThreeVector(fxcenter-fiberDiameter/2.+gapX/2.-0.2*mm,fycenter+ fiberDistance- 4*mm,fzcenter+fiberLenght3/2. -gapZ/2.+0.1*mm),
                          logicWLSGap3,//logicWLSGap1
                          "WLSfiber_quad",
                          logicMother,
                          false,
                          fibNum+Findex,
                      fCheckOverlaps);
                      
                      
        // triplet
        
        physiWLSGap_trip[i]=
        new G4PVPlacement(0,
                          G4ThreeVector(fxcenter-fiberDiameter/2.+gapX/2.-0.2*mm -OffsetForConstruction,fycenter - 4*mm + 2*mm,fzcenter+fiberLenght3/2. -gapZ/2.+0.1*mm),
                          logicWLSGap2_trip,//logicWLSGap1
                          "WLSfiber_trip",
                          logicMother,
                          false,
                          fibNum+Findex,
                      fCheckOverlaps);

        physiWLSGap_trip[i]=
        new G4PVPlacement(0,
                          G4ThreeVector(fxcenter-fiberDiameter/2.+gapX/2.-0.2*mm -OffsetForConstruction,fycenter+ fiberDistance- 4*mm + 2*mm,fzcenter+fiberLenght3/2. -gapZ/2.+0.1*mm),
                          logicWLSGap2_trip,//logicWLSGap1
                          "WLSfiber_trip",
                          logicMother,
                          false,
                          fibNum+Findex,
                      fCheckOverlaps);
      }


      if(zCalo == 2){
        
        // quadruplet
        
        physiWLSGap[i]=
        new G4PVPlacement(0,
                          G4ThreeVector(fxcenter+fiberDiameter/2.-gapX/2.+0.2*mm,fycenter,fiberLenght0/2.-gapZ+0.1*mm +3*30*mm),
                          logicWLSGap0,//logicWLSGap1
                          "WLSfiber_quad",
                          logicMother,
                          false,
                          fibNum+Findex,
                      fCheckOverlaps);

        physiWLSGap[i]=
        new G4PVPlacement(0,
                          G4ThreeVector(fxcenter+fiberDiameter/2.-gapX/2.+0.2*mm,fycenter+ fiberDistance,fiberLenght0/2.-gapZ+0.1*mm +3*30*mm),
                          logicWLSGap0,//logicWLSGap1
                          "WLSfiber_quad",
                          logicMother,
                          false,
                          fibNum+Findex,
                      fCheckOverlaps);
                      
        // triplet
        
        physiWLSGap_trip[i]=
        new G4PVPlacement(0,
                          G4ThreeVector(fxcenter+fiberDiameter/2.-gapX/2.+0.2*mm -OffsetForConstruction,fycenter,fiberLenght0/2.-gapZ+0.1*mm +3*30*mm),
                          logicWLSGap3_trip,//logicWLSGap1
                          "WLSfiber_trip",
                          logicMother,
                          false,
                          fibNum+Findex,
                      fCheckOverlaps);

        physiWLSGap_trip[i]=
        new G4PVPlacement(0,
                          G4ThreeVector(fxcenter+fiberDiameter/2.-gapX/2.+0.2*mm -OffsetForConstruction,fycenter+ fiberDistance,fiberLenght0/2.-gapZ+0.1*mm +3*30*mm),
                          logicWLSGap3_trip,//logicWLSGap1
                          "WLSfiber_trip",
                          logicMother,
                          false,
                          fibNum+Findex,
                      fCheckOverlaps);
        
        
      }
      }//for loop for fibers*/
      
      
      /*
      nbOfGaps = 1;  //5
      
      for(G4int i=0;i<nbOfGaps;i++) {
      

      G4double Gxcenter = -calorX/2.+gapX/2.+ (gapX+absX)*i + nBoxloop*totalCalorThickness;
      G4double Gycenter=gapY/2.+yCalo*gapY;
      G4double Gzcenter=gapZ/2.+zCalo*gapZ;
      G4int Gindex =zCalo*nbOfYCalo*nbOfGaps + yCalo*nbOfGaps + nBoxloop*nbOfZCalo*nbOfYCalo*nbOfGaps + i;
     
     
     */
	/*
      if(i<2){

        if (yCalo < 2){
        	// t0
          physiGap[Gindex] = new G4PVPlacement(0,
            G4ThreeVector(Gxcenter,Gycenter,Gzcenter-gapZ),
            logict0,
            "physiGap_t0",
            logicMother,
            false,
            gapNum+Gindex,
            fCheckOverlaps);

        }
        }*/
        
        
        /*
        
        if(zCalo == 0){ 
        
        	physiGap[Gindex] = new G4PVPlacement(0,
            G4ThreeVector(Gxcenter,Gycenter,Gzcenter-gapZ),
            logict0,
            "physiGap_t0",
            logicMother,
            false,
            gapNum+Gindex,
            fCheckOverlaps);
            
          // t1
          physiGap[Gindex] = new G4PVPlacement(0,
            G4ThreeVector(Gxcenter,Gycenter,Gzcenter),
            logict1,
            "physiGap_t1",
            logicMother,
            false,
            gapNum+Gindex,
            fCheckOverlaps);
          // t4
          physiGap[Gindex] = new G4PVPlacement(0,
            G4ThreeVector(Gxcenter - OffsetForConstruction,Gycenter,Gzcenter),
            logict4,
            "physiGap_t4",
            logicMother,
            false,
            gapNum+Gindex,
            fCheckOverlaps);        
            
        }

        if(zCalo == 1){
          // t2
          physiGap[Gindex] = new G4PVPlacement(0,
            G4ThreeVector(Gxcenter,Gycenter,Gzcenter),
            logict2,
            "physiGap_t2",
            logicMother,
            false,
            gapNum+Gindex,
            fCheckOverlaps);
            
          // t5
          physiGap[Gindex] = new G4PVPlacement(0,
            G4ThreeVector(Gxcenter - OffsetForConstruction,Gycenter,Gzcenter),
            logict5,
            "physiGap_t5",
            logicMother,
            false,
            gapNum+Gindex,
            fCheckOverlaps);
        }

        if(zCalo == 2){
          // t3 
            physiGap[Gindex] = new G4PVPlacement(0,
            G4ThreeVector(Gxcenter,Gycenter,Gzcenter),
            logict3,
            "physiGap_t3",
            logicMother,
            false,
            gapNum+Gindex,
            fCheckOverlaps);
            
          // t6
         
          physiGap[Gindex] = new G4PVPlacement(0,
            G4ThreeVector(Gxcenter - OffsetForConstruction,Gycenter,Gzcenter),
            logict6,
            "physiGap_t6",
            logicMother,
            false,
            gapNum+Gindex,
            fCheckOverlaps);
        }
           
      }//for loop for gaps
      
      */
      
  /* G4int Pindex0=0;
    nbOfPoli = 2; //15
    
    for(G4int i=0;i<nbOfPoli;i++) {

      G4double Pxcenter =-calorX/2.+absX+(gapX+absX)*i + nBoxloop*totalCalorThickness;
      G4double Pycenter=poliY/2.+yCalo*poliY;
      G4double Pzcenter=poliZ/2.+absZ;
      Pindex0=xCalo*nbOfXCalo*nbOfPoli + i; //nBoxloop*nbOfZCalo*nbOfYCalo*nbOfAbsor + i;

             physiPoli[Pindex0] = new G4PVPlacement(0,
             G4ThreeVector(Pxcenter - poliX,Pycenter,Pzcenter), // edited by F.B & S. M
             logicPoli,
             "physiPoli",
             logicMother,
             false,
             poliNum+Pindex0,
             fCheckOverlaps);
      //indexem0=index0;
       
    }*/
    
    
    /*
}//chiudo zCalo

}//chiudo yCalo

*/
/*
  G4int index0=0;
    nbOfAbsor = 5; //15
    for(G4int i=0;i<nbOfAbsor;i++) {

      G4double Axcenter =-calorX/2.+absX/2.+gapX+(gapX+absX)*i + nBoxloop*totalCalorThickness;
      G4double Aycenter=absY/2.;
      G4double Azcenter=absZ/2.;
      index0=xCalo*nbOfXCalo*nbOfAbsor + i; //nBoxloop*nbOfZCalo*nbOfYCalo*nbOfAbsor + i;

             physiAbs[index0] = new G4PVPlacement(0,
             G4ThreeVector(Axcenter,Aycenter,Azcenter),
             logicAbs,
             "physiAbsorber",
             logicMother,
             false,
             absNum+index0,
             fCheckOverlaps);
      //indexem0=index0;
        
    }*/
    
    
    /*
    
}//chiudo xCalo

}
*/

 /*
for(int uu=0;uu<nbOfXCalo;uu++)
  // for(int uu=0;uu<1;uu++)
   {


	physiBox1[uu] = new G4PVPlacement(0,
					  G4ThreeVector(uu*BoxSizeX1,0,0),
					     logicBox1,
					     "physiBox1",
					     logicWorld,
					     false,
					     uu,
					     fCheckOverlaps);

   }
  //}
  */



  /*
  solidBrick= new G4Box("solidBrick",12.*cm,12.*cm,12.*cm);// 24 cm thick in beam direction
  logicBrick = new G4LogicalVolume(solidBrick,brickMaterial,"logicBrick");
  physiBrick =
      new G4PVPlacement(0,
                        //G4ThreeVector(calorThickness*(nbOfXCalo+0.5),(calorSizeYZ*(nbOfYCalo-1))-calorSizeYZ/2.,(calorSizeYZ*(nbOfZCalo-1))-calorSizeYZ/2.),
			G4ThreeVector(52.5*cm,(calorSizeYZ*(nbOfYCalo-1))-calorSizeYZ/2.+8*cm,(calorSizeYZ*(nbOfZCalo-1))-calorSizeYZ/2.),
                        logicBrick,
                        "physiBrick",
                        logicWorld,
                        false,
                        10000+500,
                        fCheckOverlaps);

  // scintillator of the muon catcher
  solidMCScin= new G4Box("solidMCScin",0.5*cm,7.5*cm,7.5*cm);//check thickness
  logicMCScin = new G4LogicalVolume(solidMCScin,MCMaterial,"logicMCScin");
  physiMCScin = new G4PVPlacement(0,
				  G4ThreeVector(52.5*cm+12.*cm+1*cm,5.5*cm,0*cm),
				  logicMCScin,
				  "physiMCScin",
				  logicWorld,
				  false,
				  10000+500,
				  fCheckOverlaps);

  */
  ////SURFACE PROPERTIES///


  //ANS
  /*
  G4OpticalSurface* paintedOptSurf = new G4OpticalSurface("paintedOptSurf");
  paintedOptSurf->SetType(dielectric_dielectric);
  paintedOptSurf->SetModel(unified);
  paintedOptSurf->SetFinish(polishedfrontpainted);

  G4OpticalSurface* groundOptSurf = new G4OpticalSurface("groundOptSurf");
  groundOptSurf->SetType(dielectric_dielectric);
  groundOptSurf->SetModel(unified);
  groundOptSurf->SetFinish(ground);
  groundOptSurf->SetSigmaAlpha(0.3);
  // groundOptSurf->SetFinish(polished);

  const G4int numPhotonEnergy=30;

  G4double photonEnergy[numPhotonEnergy]={2.066*eV,2.101*eV,2.137*eV,2.175*eV,2.214*eV,
                                          2.254*eV,2.296*eV,2.339*eV,2.384*eV,2.431*eV,
                                          2.455*eV,2.479*eV,2.530*eV,2.556*eV,2.583*eV,
                                          2.610*eV,2.637*eV,2.695*eV,2.724*eV,2.755*eV,
                                          2.786*eV,2.818*eV,2.850*eV,2.883*eV,2.917*eV,
                                          2.952*eV,2.987*eV,3.024*eV,3.099*eV,3.179*eV};

  G4double reflectivity[numPhotonEnergy] = {0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,
                                            0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,
                                            0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9};
  // G4double reflectivity[numPhotonEnergy] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
  //                                           0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
  //                                           0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};

  G4MaterialPropertiesTable* paintedOptSurfTable = new G4MaterialPropertiesTable();
  paintedOptSurfTable->AddProperty("REFLECTIVITY",photonEnergy,reflectivity,numPhotonEnergy);
  paintedOptSurf->SetMaterialPropertiesTable(paintedOptSurfTable);

  G4double specularLobe[numPhotonEnergy] = {0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,
                                            0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,
                                            0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3};
  G4double specularSpike[numPhotonEnergy] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
                                             0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
                                             0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
  G4double backScatter[numPhotonEnergy] = {0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,
                                           0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,
                                           0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3};

  G4MaterialPropertiesTable* groundOptSurfTable = new G4MaterialPropertiesTable();
  groundOptSurfTable->AddProperty("SPECULARLOBECONSTANT",photonEnergy,specularLobe,numPhotonEnergy);
  groundOptSurfTable->AddProperty("SPECULARSPIKECONSTANT",photonEnergy,specularSpike,numPhotonEnergy);
  groundOptSurfTable->AddProperty("BACKSCATTERCONSTANT",photonEnergy,backScatter,numPhotonEnergy);
  groundOptSurf->SetMaterialPropertiesTable(groundOptSurfTable);




  G4LogicalBorderSurface* scintiAbsorberSurface[maxTiles*2];
  G4LogicalBorderSurface* scintiWorldSurface[maxTiles];

  for (G4int i=0;i<nbOfGaps*nbOfXCalo*nbOfYCalo*nbOfZCalo;i++) {
    scintiWorldSurface[i]= new G4LogicalBorderSurface("scintiWorldSurface",physiGap[i],physiWorld,paintedOptSurf);

    for(G4int j=0;j<2;j++)
      scintiAbsorberSurface[2*i+j] = new G4LogicalBorderSurface("scintiAbsorberSurface",physiGap[i],physiAbs[i+j],paintedOptSurf);
  }

  G4LogicalBorderSurface* scintiHoleSurface[maxTiles*maxFibers];
  G4LogicalBorderSurface* holeScintiSurface[maxTiles*maxFibers];
  for (G4int i=0;i<nbOfGaps*nbOfXCalo*nbOfYCalo*nbOfZCalo;i++) {
    for (G4int j=0;j<nbOfFibers;j++) {
      scintiHoleSurface[nbOfFibers*i+j] = new G4LogicalBorderSurface("scintiHoleSurface",physiGap[i],physiWLSGap[j],groundOptSurf);
      holeScintiSurface[nbOfFibers*i+j] = new G4LogicalBorderSurface("holeScintiSurface",physiWLSGap[j],physiGap[i],groundOptSurf);
    }
  }


  G4LogicalBorderSurface* scintiSensorPlaneSurface[maxTiles];

  for (G4int i=0;i<nbOfXCalo*nbOfYCalo*nbOfZCalo;i++) {
    scintiSensorPlaneSurface[i] =  new G4LogicalBorderSurface("scintiSensorPlaneSurface",physiGap[nbOfGaps*(i+1)-1],physiSensorPlane[i],paintedOptSurf);
  }
  //-------------------------------------------------------------------//
  */
  
  /*
  PrintCalorParameters();
}

*/

G4LogicalVolume * DetectorConstruction::ConstructWLSFiber(G4double length) {

 solidFiber = new G4Tubs("solidCore",0.0,(fiberDiameter-fiberDiameter*6/100.)/2.,
			  length/2.,0.0*deg,360.*deg);

  solidClad1 = new G4Tubs("solidCladding1",0.0,(fiberDiameter-fiberDiameter*3/100.)/2.,
			  length/2.,0.0*deg,360.*deg);

  solidClad2 = new G4Tubs("solidCladding2",0.0,fiberDiameter/2.,
			  length/2.,0.0*deg,360.*deg);

  solidHole = new G4Tubs("solidHole",0.0,(fiberDiameter+0.1*mm)/2.,
                         length/2.,0.0*deg,360.*deg);

  logicFiber = new G4LogicalVolume(solidFiber,fiberMaterial,"logicCore");
  logicClad1 = new G4LogicalVolume(solidClad1,cladding1Material,"logicCladding1");
  logicClad2 = new G4LogicalVolume(solidClad2,cladding2Material,"logicCladding2");
  logicHole = new G4LogicalVolume(solidHole,worldMaterial,"logicHole");

  physiFiber = new G4PVPlacement(0,G4ThreeVector(0,0,0),logicFiber,"physiCore",
                                 logicClad1,false,1,fCheckOverlaps);

  physiClad1 = new G4PVPlacement(0,G4ThreeVector(0,0,0),logicClad1,"physiCladding1",
                                 logicClad2,false,1,fCheckOverlaps);

  physiClad2 = new G4PVPlacement(0,G4ThreeVector(0,0,0),logicClad2,"physiCladding2",
                                 logicHole,false,1,fCheckOverlaps);



  return logicHole;
}

void DetectorConstruction::PrintCalorParameters()
{

  G4cout << "\n-------------------------------------------------------------"
	 << "\n ---> The calorimeter is " << nbOfAbsor << " tiles of:";

  G4cout << "\n \t" << std::setw(12) << absorMaterial->GetName() <<": "
	 << std::setw(6) << G4BestUnit(absX,"Length");

  G4cout << "\n-----------------------------------------------------------\n";
  G4cout << "\n" << absorMaterial << G4endl;

  G4cout << "\n-----------------------------------------------------------\n";

  G4cout << "\n-------------------------------------------------------------"
	 << "\n ---> The calorimeter is " << nbOfGaps << " tiles of:";

  G4cout << "\n \t" << std::setw(12) << gapMaterial->GetName() <<": "
	 << std::setw(6) << G4BestUnit(gapZ,"Length");

  G4cout << "\n-----------------------------------------------------------\n";
  G4cout << "\n" << gapMaterial << G4endl;

  G4cout << "\n-----------------------------------------------------------\n";


  G4cout << "\n-------------------------------------------------------------"
	 << "\n ---> For a total of " << nbOfPlanes << " tiles."<< G4endl;
  G4cout << "\n-----------------------------------------------------------\n";

}

void DetectorConstruction::SetAbsorThickness(G4double val) {

  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetAbsorThickness: thickness "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  absX = val;
}

void DetectorConstruction::SetGapThickness(G4double val) {

  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetAbsorThickness: thickness "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  gapX = val;
}

void DetectorConstruction::SetCalorSizeYZ(G4double val) {

  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetCalorSizeYZ: thickness "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  absY = val;
}

void DetectorConstruction::ConstructSDandField() {

  SipmSensitiveDetector* sipmSD
    = new SipmSensitiveDetector("SipmSD");
    SetSensitiveDetector("logicSipm",sipmSD);

  GapSensitiveDetector* gapSD
    = new GapSensitiveDetector("GapSD");  
  G4SDManager::GetSDMpointer()->AddNewDetector(gapSD);
   
  for (int j = 0; j < N_tile_t0; j++) {	
  	SetSensitiveDetector(Form("logict0[%d]",j),gapSD);
  }
  for (int j = 0; j < N_tile_LCM; j++) {	
  	SetSensitiveDetector(Form("logict1[%d]",j),gapSD);
  	SetSensitiveDetector(Form("logict2[%d]",j),gapSD);
  	SetSensitiveDetector(Form("logict3[%d]",j),gapSD);
  	SetSensitiveDetector(Form("logict4[%d]",j),gapSD);
  	SetSensitiveDetector(Form("logict5[%d]",j),gapSD);
  	SetSensitiveDetector(Form("logict6[%d]",j),gapSD);
  }

  AbsSensitiveDetector* absSD
    = new AbsSensitiveDetector("AbsSD");
  SetSensitiveDetector("logicAbsorber",absSD);
  
  MuonSensitiveDetector* muonSD
    = new MuonSensitiveDetector("muonSD");
  SetSensitiveDetector("logicMuon1", muonSD);
  SetSensitiveDetector("logicMuon2", muonSD);
  


  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  G4AutoDelete::Register(fMagFieldMessenger);
}


// ************ EXPERIMENTAL **********








void DetectorConstruction::ConstructTiles() {

  G4VisAttributes* poliAttr= new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  G4VisAttributes* motherAttr= new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  G4VisAttributes* gapAttr= new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  
  
  G4VisAttributes* absAttr= new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  G4VisAttributes* fiberAttr= new G4VisAttributes(G4Colour(0.0,0.8,0.0));
  G4VisAttributes* clad1Attr= new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  G4VisAttributes* clad2Attr= new G4VisAttributes(G4Colour(0.0,0.5,0.0));
  G4VisAttributes* sipmAttr= new G4VisAttributes(G4Colour(0.8,0.8,0.8));
  G4VisAttributes* sPlaneAttr= new G4VisAttributes(G4Colour(0.1,0.6,0.1));

  clad1Attr->SetVisibility(true);
  clad1Attr->SetForceSolid(true);

  clad2Attr->SetVisibility(true);
  clad2Attr->SetForceSolid(true);

  absAttr->SetForceSolid(true);
  gapAttr->SetVisibility(true);
  poliAttr->SetVisibility(true);

  motherAttr->SetVisibility(false);

  fiberAttr->SetVisibility(true);
  fiberAttr->SetForceSolid(true);

  sipmAttr->SetForceSolid(true);

  
  // ---------- Solids ----------
  
  // ----- Randomization of tile thicknesses -----
  
  
  //CLHEP::RandGauss rnd = CLHEP::RandGauss( *(CLHEP::HepRandom::getTheEngine()) , 0.0 , 1.0 );  
  
  
  
  
  for (int j = 0; j < N_tile_t0; j++) {

  	tile_thickness_eff = RandomThickness(tile_thickness_nom);
  	G4cout << "thick eff = " << tile_thickness_eff << G4endl;
	t0_thick.push_back(tile_thickness_eff);
  	THICK_OFFSET = tile_thickness_eff - tile_thickness_nom;

  	G4Trd* tile0 = new G4Trd("tile0", tile_thickness_eff/2., tile_thickness_eff/2., t0_min_base/2., t0_maj_base/2., h/2.);
  	G4Box* groove = new G4Box("innerfiber", grooveDiameter/2. + THICK_OFFSET/4., grooveDiameter/2., h/2.); // asymmetrical groove !! THICKNESS OFFSET
  	G4SubtractionSolid* hollowGroove_t0 = new G4SubtractionSolid("hollow groove of tile0", tile0, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  	
	//t0[j] = new G4SubtractionSolid(Form("t0[%d]",j), hollowGroove_t0, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
	t0.push_back( new G4SubtractionSolid(Form("t0[%d]",j), hollowGroove_t0, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.)) );
	
	
	/*
	delete tile0;
	delete groove;
	delete hollowGroove_t0;*/
	 	
  }
	
  //Filling vector of tiles thickness after every creation of N tile of the same type
  random_thickness.push_back(t0_thick);
  
  
  // tiles in the calorimeters
  // quadruplets
  // t1
  for (int j = 0; j < N_tile_LCM; j++) {
  
  	tile_thickness_eff = RandomThickness(tile_thickness_nom);
  	G4cout << "thick eff = " << tile_thickness_eff << G4endl;
	t1_thick.push_back(tile_thickness_eff);
  	THICK_OFFSET = tile_thickness_eff - tile_thickness_nom;
  	
  	G4Trd* tile1 = new G4Trd("tile1", tile_thickness_eff/2., tile_thickness_eff/2., t1_min_base/2., t1_maj_base/2., h/2.);
  	
  	G4Box* groove = new G4Box("innerfiber", grooveDiameter/2. + THICK_OFFSET/4., grooveDiameter/2., h/2.);
  	G4SubtractionSolid *hollowGroove_t1_1 = new G4SubtractionSolid("hollow groove of tile1 - 1", tile1, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  	G4SubtractionSolid *hollowGroove_t1_2 = new G4SubtractionSolid("hollow groove of tile1 - 2", hollowGroove_t1_1, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  	G4SubtractionSolid *hollowGroove_t1_3 = new G4SubtractionSolid("hollow groove of tile1 - 3", hollowGroove_t1_2, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2)/2., -fiberDistance/2.,0.));
  	t1.push_back( new G4SubtractionSolid(Form("t1[%d]",j), hollowGroove_t1_3, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.)) );
  	/*
  	delete tile1;
	delete groove;
	delete hollowGroove_t1_1;
	delete hollowGroove_t1_2;
	delete hollowGroove_t1_3;*/
  	
  }
	random_thickness.push_back(t1_thick);
  // t2
  for (int j = 0; j < N_tile_LCM; j++) {
  	
  	tile_thickness_eff = RandomThickness(tile_thickness_nom);
  	G4cout << "thick eff = " << tile_thickness_eff << G4endl;
	t2_thick.push_back(tile_thickness_eff);
  	THICK_OFFSET = tile_thickness_eff - tile_thickness_nom;
  	
  	G4Trd* tile2 = new G4Trd("tile2", tile_thickness_eff/2., tile_thickness_eff/2., t2_min_base/2., t2_maj_base/2., h/2.);
  	G4Box* groove = new G4Box("innerfiber", grooveDiameter/2. + THICK_OFFSET/4., grooveDiameter/2., h/2.); // asymmetrical groove !! THICKNESS OFFSET
  	G4SubtractionSolid *hollowGroove_t2_1 = new G4SubtractionSolid("hollow groove of tile2 - 1", tile2, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  	G4SubtractionSolid *hollowGroove_t2_2 = new G4SubtractionSolid("hollow groove of tile2 - 2", hollowGroove_t2_1, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  	G4SubtractionSolid *hollowGroove_t2_3 = new G4SubtractionSolid("hollow groove of tile2 - 3", hollowGroove_t2_2, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  	G4SubtractionSolid *hollowGroove_t2_4 = new G4SubtractionSolid("hollow groove of tile2 - 4", hollowGroove_t2_3, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  	G4SubtractionSolid * hollowGroove_t2_5 = new G4SubtractionSolid("hollow groove of tile2 - 5", hollowGroove_t2_4, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. -2*mm,0.));
  	t2.push_back( new G4SubtractionSolid(Form("t2[%d]",j), hollowGroove_t2_5, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance -2*mm,0.)) );
  	/*
  	delete tile2;
	delete groove;
	delete hollowGroove_t2_1;
	delete hollowGroove_t2_2;
	delete hollowGroove_t2_3;
	delete hollowGroove_t2_4;
	delete hollowGroove_t2_5;*/
  	
  }
	random_thickness.push_back(t2_thick);
  // t3
  for (int j = 0; j < N_tile_LCM; j++) {
  
  	tile_thickness_eff = RandomThickness(tile_thickness_nom);
  	G4cout << "thick eff = " << tile_thickness_eff << G4endl;
	t3_thick.push_back(tile_thickness_eff);
  	THICK_OFFSET = tile_thickness_eff - tile_thickness_nom;
  	
  	G4Trd* tile3 = new G4Trd("tile3", tile_thickness_eff/2., tile_thickness_eff/2., t3_min_base/2., t3_maj_base/2., h/2.);
  	G4Box* groove = new G4Box("innerfiber", grooveDiameter/2. + THICK_OFFSET/4., grooveDiameter/2., h/2.); // asymmetrical groove !! THICKNESS OFFSET
  	
  	G4SubtractionSolid *hollowGroove_t3_1 = new G4SubtractionSolid("hollow groove of tile3 - 1", tile3, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  	G4SubtractionSolid *hollowGroove_t3_2 = new G4SubtractionSolid("hollow groove of tile3 - 2", hollowGroove_t3_1, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  	G4SubtractionSolid *hollowGroove_t3_3 = new G4SubtractionSolid("hollow groove of tile3 - 3", hollowGroove_t3_2, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  	G4SubtractionSolid *hollowGroove_t3_4 = new G4SubtractionSolid("hollow groove of tile3 - 4", hollowGroove_t3_3, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  	G4SubtractionSolid * hollowGroove_t3_5 = new G4SubtractionSolid("hollow groove of tile3 - 5", hollowGroove_t3_4, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. -2*mm,0.));
  	G4SubtractionSolid * hollowGroove_t3_6 = new G4SubtractionSolid("hollow groove of tile3 - 6", hollowGroove_t3_5, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance -2*mm,0.));
  	G4SubtractionSolid * hollowGroove_t3_7 = new G4SubtractionSolid("hollow groove of tile3 - 7", hollowGroove_t3_6, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + 2*mm,0.));
  	t3.push_back( new G4SubtractionSolid(Form("t3[%d]",j), hollowGroove_t3_7, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance + 2*mm,0.)) );
  	/*
  	delete tile3;
	delete groove;
	delete hollowGroove_t3_1;
	delete hollowGroove_t3_2;
	delete hollowGroove_t3_3;
	delete hollowGroove_t3_4;
	delete hollowGroove_t3_5;
	delete hollowGroove_t3_6;
	delete hollowGroove_t3_7;*/
  }
	random_thickness.push_back(t3_thick);
  // triplets
  
  // t4
  for (int j = 0; j < N_tile_LCM; j++) {
  	
  	tile_thickness_eff = RandomThickness(tile_thickness_nom);
  	G4cout << "thick eff = " << tile_thickness_eff << G4endl;
    t4_thick.push_back(tile_thickness_eff);
  	THICK_OFFSET = tile_thickness_eff - tile_thickness_nom;
  		
  	G4Trd* tile4 = new G4Trd("tile4", tile_thickness_eff/2., tile_thickness_eff/2., t4_min_base/2., t4_maj_base/2., h/2.);
  	G4Box* groove = new G4Box("innerfiber", grooveDiameter/2. + THICK_OFFSET/4., grooveDiameter/2., h/2.); // asymmetrical groove !! THICKNESS OFFSET
  	
  	G4SubtractionSolid *hollowGroove_t4 = new G4SubtractionSolid(Form("hollow groove of tile4 [%d]",j), tile4, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  	t4.push_back( new G4SubtractionSolid(Form("t4[%d]",j), hollowGroove_t4, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.)) );
  	/*
  	delete tile4;
	delete groove;
	delete hollowGroove_t4; */
  }
	random_thickness_trip.push_back(t4_thick);
  
  // t5
  for (int j = 0; j < N_tile_LCM; j++) {
  	
  	tile_thickness_eff = RandomThickness(tile_thickness_nom);
  	G4cout << "thick eff = " << tile_thickness_eff << G4endl;
	t5_thick.push_back(tile_thickness_eff);
  	THICK_OFFSET = tile_thickness_eff - tile_thickness_nom;
  	
  	G4Trd* tile5 = new G4Trd("tile5", tile_thickness_eff/2., tile_thickness_eff/2., t5_min_base/2., t5_maj_base/2., h/2.);
  	G4Box* groove = new G4Box("innerfiber", grooveDiameter/2. + THICK_OFFSET/4., grooveDiameter/2., h/2.); // asymmetrical groove !! THICKNESS OFFSET
  	
  	G4SubtractionSolid *hollowGroove_t5_1 = new G4SubtractionSolid("hollow groove of tile4 - 1", tile5, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  	G4SubtractionSolid *hollowGroove_t5_2 = new G4SubtractionSolid("hollow groove of tile4 -2", hollowGroove_t5_1, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  	G4SubtractionSolid *hollowGroove_t5_3 = new G4SubtractionSolid("hollow groove of tile4 -3", hollowGroove_t5_2, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  	t5.push_back( new G4SubtractionSolid(Form("t5[%d]",j), hollowGroove_t5_3, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.)) );
  	/*
  	delete tile5;
	delete groove;
	delete hollowGroove_t5_1;
	delete hollowGroove_t5_2;
	delete hollowGroove_t5_3;*/
	}
	random_thickness_trip.push_back(t5_thick);
	
  // t6
  for (int j = 0; j < N_tile_LCM; j++) {
  	
  	tile_thickness_eff = RandomThickness(tile_thickness_nom);
  	G4cout << "thick eff = " << tile_thickness_eff << G4endl;
	t6_thick.push_back(tile_thickness_eff);
  	THICK_OFFSET = tile_thickness_eff - tile_thickness_nom;
  	
  	G4Trd* tile6 = new G4Trd("tile6", tile_thickness_eff/2., tile_thickness_eff/2., t6_min_base/2., t6_maj_base/2., h/2.);
  	G4Box* groove = new G4Box("innerfiber", grooveDiameter/2. + THICK_OFFSET/4., grooveDiameter/2., h/2.); // asymmetrical groove !! THICKNESS OFFSET
  	
  	G4SubtractionSolid *hollowGroove_t6_1 = new G4SubtractionSolid("hollow groove of tile6 - 1", tile6, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  	G4SubtractionSolid *hollowGroove_t6_2 = new G4SubtractionSolid("hollow groove of tile6 -2", hollowGroove_t6_1, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  	G4SubtractionSolid *hollowGroove_t6_3 = new G4SubtractionSolid("hollow groove of tile6 -3", hollowGroove_t6_2, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2.,0.));
  	G4SubtractionSolid *hollowGroove_t6_4 = new G4SubtractionSolid("hollow groove of tile6 -4", hollowGroove_t6_3, groove, 0, G4ThreeVector(-(tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance,0.));
  	G4SubtractionSolid *hollowGroove_t6_5 = new G4SubtractionSolid("hollow groove of tile6 -5", hollowGroove_t6_4, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. + fiberDistance - 2*mm,0.));
  
  	t6.push_back( new G4SubtractionSolid(Form("t6[%d]",j), hollowGroove_t6_5, groove, 0, G4ThreeVector((tile_thickness_nom - grooveDiameter + THICK_OFFSET/2.)/2., -fiberDistance/2. - 2*mm,0.)) );
  	
  	// some problems happens if you un-comment
  	/*
  	delete tile6;
	delete groove;
	delete hollowGroove_t6_1;
	delete hollowGroove_t6_2;
	delete hollowGroove_t6_3;
	delete hollowGroove_t6_4;
	delete hollowGroove_t6_5;
	*/
  }
	
	random_thickness_trip.push_back(t6_thick);
  
  
  /*
//polietilene borato
  G4Box *outerPoli = new G4Box("tile 0",poliX/2., poliY/2., poliZ/2.);
  G4Box* groove_Poli = new G4Box("outer groove",grooveX/2., grooveY/2., grooveZ/2.); 
  //G4Tubs *innerfiber2 = new G4Tubs ("innerfiber2", 0, fiberDiameter/2.+0.5*mm, poliZ/2., 0, 2*M_PI*rad);
  G4SubtractionSolid *hollowPoli = new G4SubtractionSolid("Hollow Poli", outerPoli, groove_Poli, 0, G4ThreeVector(-(poliX/2.-0.55*mm),0.,0.));
  G4SubtractionSolid *poliFinal = new G4SubtractionSolid("Hollow Poli", hollowPoli, groove_Poli, 0, G4ThreeVector(+(poliX/2.-0.55*mm),0.,0.));*/
  
  /*
  // FIX ME: IT was like this before
  //polietilene borato
  G4Box *outerPoli = new G4Box("tile 0",poliX/2., poliY/2., poliZ/2.);
  G4Box* grooveBPE = new G4Box("outer groove",grooveX/2., grooveY/2., grooveZ/2.); 
  //G4Tubs *innerfiber2 = new G4Tubs ("innerfiber2", 0, fiberDiameter/2.+0.5*mm, poliZ/2., 0, 2*M_PI*rad);
   G4SubtractionSolid *hollowPoli = new G4SubtractionSolid("Hollow Poli", outerPoli, grooveBPE, 0, G4ThreeVector(-(poliX/2.-0.55*mm),0.,0.));
  G4SubtractionSolid *poliFinal = new G4SubtractionSolid("Hollow Poli", hollowPoli, grooveBPE, 0, G4ThreeVector(+(poliX/2.-0.55*mm),0.,0.));*/
  
  
  // 
 
  solidAbs = new G4Box("solidAbsorber",absX/2.,absY/2.,absZ/2.);
  solidGap = new G4Box("solidGap",gapX/2.,gapY/2.,gapZ/2.);
  solidPoli = new G4Box("solidPoli", poliX/2., poliY/2., poliZ/2.); 
  //scatola esterna
  solidMother = new G4Box("solidMother",motherX/2.,motherY/2.,motherZ/2.);

  G4double sipmDimX=0.5*mm;
  G4double sipmDimY=fiberDiameter/2.;
  G4double sipmDimZ=fiberDiameter/2.;
  G4double airCushionDimX = sipmDimX/10.;
  
  
  
  



	
	
  
  
  
  

  solidSensorPlane = new G4Box("solidSensorPlane",sensPlaneThickness/2.,absY/2.,absZ/2.);
  solidSipm = new G4Box("solidSipm",sipmDimX,sipmDimY,sipmDimZ);
  solidAirCushion = new G4Box("solidAirCushion",airCushionDimX,sipmDimY,sipmDimZ);
  
  
  
  // ---------- Logical Volumes ----------
  
  logicMother = new G4LogicalVolume(solidMother,worldMaterial,"logicMother");
  logicAbs = new G4LogicalVolume(solidAbs,absorMaterial,"logicAbsorber");
    
  // tiles in t0 layer
  for (int j = 0; j < N_tile_t0; j++) {
  	logict0.push_back( new G4LogicalVolume(t0[j],gapMaterial,Form("logict0[%d]",j)) );
  }
  // tiles in the calorimeters
  for (int j = 0; j < N_tile_LCM; j++) {
  	logict1.push_back( new G4LogicalVolume(t1[j],gapMaterial,Form("logict1[%d]",j)) );
  	logict2.push_back( new G4LogicalVolume(t2[j],gapMaterial,Form("logict2[%d]",j)) );
  	logict3.push_back( new G4LogicalVolume(t3[j],gapMaterial,Form("logict3[%d]",j)) );
  	logict4.push_back( new G4LogicalVolume(t4[j],gapMaterial,Form("logict4[%d]",j)) );
  	logict5.push_back( new G4LogicalVolume(t5[j],gapMaterial,Form("logict5[%d]",j)) );
  	logict6.push_back( new G4LogicalVolume(t6[j],gapMaterial,Form("logict6[%d]",j)) );
  }
  // logicPoli = new G4LogicalVolume(poliFinal, poliMaterial, "poliFinal");


  //logicWLSAbs = ConstructWLSFiber(calorSizeYZ);
  // quadruplet
  logicWLSGap0 = ConstructWLSFiber(fiberLenght0);
  logicWLSGap1 = ConstructWLSFiber(fiberLenght1);
  logicWLSGap2 = ConstructWLSFiber(fiberLenght2);
  logicWLSGap3 = ConstructWLSFiber(fiberLenght3);
  // triplet
  logicWLSGap1_trip = ConstructWLSFiber(fiberLenght1);
  logicWLSGap2_trip = ConstructWLSFiber(fiberLenght2);
  logicWLSGap3_trip = ConstructWLSFiber(fiberLenght3); // CHANGE THIS NAME


  
	
  
  
  logicSensorPlane = new G4LogicalVolume(solidSensorPlane,sensorPlaneMaterial,"logicSensorPlane");
  logicSipm = new G4LogicalVolume(solidSipm,sipmMaterial,"logicSipm");
  logicAirCushion = new G4LogicalVolume(solidAirCushion,airCushionMaterial,"logicAirCushion");





  // Set Vis Attributes
  logicAbs->SetVisAttributes(absAttr); 
  for (int j = 0; j < N_tile_t0; j++) {
  	logict0[j]->SetVisAttributes(gapAttr);
  }
  for (int j = 0; j < N_tile_LCM; j++) {
  	logict1[j]->SetVisAttributes(gapAttr);
  	logict2[j]->SetVisAttributes(gapAttr);
  	logict3[j]->SetVisAttributes(gapAttr);
  	logict4[j]->SetVisAttributes(gapAttr);
  	logict5[j]->SetVisAttributes(gapAttr);
  	logict6[j]->SetVisAttributes(gapAttr);
  }
  
  
  
  //logicPoli->SetVisAttributes(poliAttr);
  logicMother->SetVisAttributes(motherAttr);

  logicWLSGap0->SetVisAttributes(fiberAttr);
  logicWLSGap1->SetVisAttributes(fiberAttr); 
  logicWLSGap2->SetVisAttributes(fiberAttr); 
  logicWLSGap3->SetVisAttributes(fiberAttr);
  // triplet
  logicWLSGap1_trip->SetVisAttributes(fiberAttr);
  logicWLSGap2_trip->SetVisAttributes(fiberAttr);
  logicWLSGap3_trip->SetVisAttributes(fiberAttr);
  
  logicFiber->SetVisAttributes(fiberAttr);
  logicClad1->SetVisAttributes(clad1Attr);
  logicClad2->SetVisAttributes(clad2Attr);
  logicSipm->SetVisAttributes(sipmAttr);
  logicSensorPlane->SetVisAttributes(sPlaneAttr);
  logicAirCushion->SetVisAttributes(clad1Attr);

} 



G4double DetectorConstruction::RandomThickness(G4double mean) {

	G4double randThick = 0*mm;
	do {
		randThick = rnd.fire(mean, sigma);		 
	} while ( randThick >= layers_dist );
	
	return randThick;
}


// ************* Edited by S. Marangoni **************************

//To generate fiber holes in one module of BPE, using its logical volume and maximum tile thickness of the corresponding quadruplet/triplet. It has to be physically positioned once for every Phi_layer.
/* void DetectorConstruction::BPEModule_holes (G4double x_thickness) {
	G4Box* straw = new G4Box("innerHole", grooveDiameter/2., grooveDiameter/2., fiberLenght0/2.);
	//BPE module to make holes for fibers
	G4Tubs* module = new G4Tubs("module",  R_t0+h+D_Ra+D_Rg, (R_t0+h+D_Ra+D_Rg) + D_Rb, layers_dist/2., (pi-th_tile)/2., th_tile );
	//Rotation to adjust position of G4box w.r.t. position of G4Tubs
	G4RotationMatrix* r_hole = new G4RotationMatrix();
	r_hole->rotateX(90*deg);
	
//	G4cout << "LINE = "<<__LINE__ << G4endl;
	// Eight holes to let all fibers go straightforward through BPE
		//t0 L
	G4SubtractionSolid* mod_hole1 = new G4SubtractionSolid ("hollow1", module, straw, r_hole, G4ThreeVector( -fiberDistance/2., (R_t0+h+D_Ra+D_Rg) + D_Rb/2., -(layers_dist+grooveDiameter+abs(x_thickness-tile_thickness_nom))/2.+x_thickness ) );
		//t0 R
	G4SubtractionSolid* mod_hole2 = new G4SubtractionSolid ("hollow2", mod_hole1, straw, r_hole, G4ThreeVector( +fiberDistance/2., (R_t0+h+D_Ra+D_Rg) + D_Rb/2., -(layers_dist+grooveDiameter+abs(x_thickness-tile_thickness_nom))/2.+x_thickness ) );
		//t1 L
	G4SubtractionSolid* mod_hole3 = new G4SubtractionSolid ("hollow3", mod_hole2, straw, r_hole, G4ThreeVector( -fiberDistance/2., (R_t0+h+D_Ra+D_Rg) + D_Rb/2., (-layers_dist+grooveDiameter+abs(x_thickness-tile_thickness_nom))/2.  ) );
		//t1 R
	G4SubtractionSolid* mod_hole4 = new G4SubtractionSolid ("hollow4", mod_hole3, straw, r_hole, G4ThreeVector( +fiberDistance/2., (R_t0+h+D_Ra+D_Rg) + D_Rb/2., (-layers_dist+grooveDiameter+abs(x_thickness-tile_thickness_nom))/2.  ) );
		//t2 L
	G4SubtractionSolid* mod_hole5 = new G4SubtractionSolid ("hollow5", mod_hole4, straw, r_hole, G4ThreeVector( -fiberDistance/2.-2*mm, (R_t0+h+D_Ra+D_Rg) + D_Rb/2., -(layers_dist+grooveDiameter+abs(x_thickness-tile_thickness_nom))/2.+x_thickness ) );
		//t2 R
	G4SubtractionSolid* mod_hole6 = new G4SubtractionSolid ("hollow6", mod_hole5, straw, r_hole, G4ThreeVector( +fiberDistance/2.-2*mm, (R_t0+h+D_Ra+D_Rg) + D_Rb/2., -(layers_dist+grooveDiameter+abs(x_thickness-tile_thickness_nom))/2.+x_thickness ) );
		//t3 L
	G4SubtractionSolid* mod_hole7 = new G4SubtractionSolid ("hollow7", mod_hole6, straw, r_hole, G4ThreeVector( -fiberDistance/2.+2*mm, (R_t0+h+D_Ra+D_Rg) + D_Rb/2., (-layers_dist+grooveDiameter+abs(x_thickness-tile_thickness_nom))/2.  ) );
		//t3 R and final solid configuration
	solidBorPoly_groove = new G4SubtractionSolid ("solidModBPE_hollows8", mod_hole7, straw, r_hole, G4ThreeVector( +fiberDistance/2.+2*mm, (R_t0+h+D_Ra+D_Rg) + D_Rb/2., (-layers_dist+grooveDiameter+abs(x_thickness-tile_thickness_nom))/2.  ) );
}*/

//To generate fiber holes (two big grooves to contain 4 fiber each) in one module of borated polyethylene, using its logical volume. To be replicated to cover a full arc

void DetectorConstruction::BPEModule() {
	//BPE module to make holes for fibers
	G4Tubs* module = new G4Tubs("module",  R_t0+h+D_Ra+D_Rg, (R_t0+h+D_Ra+D_Rg) + D_Rb, layers_dist/2., (pi-th_tile)/2., th_tile );
	//Creating G4Box that has to be used to dig channels inside BPE module, which is a G4Tubs
	G4Box* channel = new G4Box("innerHole", grooveDiameter/2.+2.*mm, layers_dist/2.+1*mm, fiberLenght0/2.);
	//Rotation to adjust position of G4Box w.r.t. position of G4Tubs
	G4RotationMatrix* r_channel = new G4RotationMatrix();
	r_channel->rotateX(90*deg);
	
	//Two channels to let 4 fibers each go straightforward through BPE
	G4SubtractionSolid* mod_channel1 = new G4SubtractionSolid( "channel1", module, channel, r_channel, G4ThreeVector( -fiberDistance/2., (R_t0+h+D_Ra+D_Rg) + D_Rb/2. , 0.) );
	
	solidBorPoly_groove = new G4SubtractionSolid ("solidBPEModule", mod_channel1, channel, r_channel, G4ThreeVector(+fiberDistance/2. , (R_t0+h+D_Ra+D_Rg) + D_Rb/2. , 0. ) );
}
 
void DetectorConstruction::ConstructIronAbs_BorPoly() {
    
    ComputeDemonstratorParameters();
    
    solidIronAbs = new G4Tubs("SolidIronAbsorber", R_t0+h, R_t0+h+D_Ra, zAbs/2., start_phi, arc_angle);
    logicIronAbs = new G4LogicalVolume(solidIronAbs, ironAbsMaterial, "LogicIronAbsorber");
    
    G4VisAttributes* ironAbsColour = new G4VisAttributes(G4Colour::Gray());
    logicIronAbs->SetVisAttributes (ironAbsColour);
    
    solidBorPoly = new G4Tubs("SolidBoratedPolyethylene", R_t0+h+D_Ra+D_Rg, (R_t0+h+D_Ra+D_Rg) + D_Rb, zAbs/2., start_phi, arc_angle);
    logicBorPoly = new G4LogicalVolume(solidBorPoly, borpolyMaterial, "LogicBoratedPolyethylene");
	
	//In order to cover one layer with a full arc of BPE, two smaller arcs on the left and on the right of the modules with fiber holes are needed.
	solidArcL = new G4Tubs("solidArc",  R_t0+h+D_Ra+D_Rg, (R_t0+h+D_Ra+D_Rg) + D_Rb, layers_dist/2., pi/4., pi/4.-5*th_tile );
	solidArcR = new G4Tubs("solidArc",  R_t0+h+D_Ra+D_Rg, (R_t0+h+D_Ra+D_Rg) + D_Rb, layers_dist/2., pi/2.+5*th_tile, pi/4.-5*th_tile );
	logicArcL = new G4LogicalVolume(solidArcL, borpolyMaterial, "logicArc");
	logicArcR = new G4LogicalVolume(solidArcR, borpolyMaterial, "logicArc");
	
	//Arc mother for replicas of modules, which is positioned between ArcL and ArcR to cover a full arc of BPE
//	solidMotherArc = new G4Tubs("solidMotherArc", R_t0+h+D_Ra+D_Rg, (R_t0+h+D_Ra+D_Rg) + D_Rb, layers_dist/2., pi/2.-5*th_tile, 10*th_tile );
//	logicMotherArc = new G4LogicalVolume(solidMotherArc, worldMaterial, "logicMotherArc");
    
    G4VisAttributes* borpolyColour = new G4VisAttributes(G4Colour::White());
    logicBorPoly->SetVisAttributes (borpolyColour);
	logicArcL->SetVisAttributes( borpolyColour );
	logicArcR->SetVisAttributes( borpolyColour );
//	logicMotherArc->SetVisAttributes(G4VisAttributes::GetInvisible());
	
	//**** Creating replicas of a single BPE module (inside MotherArc LV) with two big holes for every Phi layer ****
//	DetectorConstruction::BPEModule();

	//Logical and physical volumes of a single replica of BPE module with fiber channels
//	logicBorPoly_groove = new G4LogicalVolume(solidBorPoly_groove, borpolyMaterial, "LogicBoratedPolyethylene");
//	logicBorPoly_groove->SetVisAttributes(borpolyColour);
//  physiBorPoly_groove = new G4PVReplica ( "PhysisBorPoly_Module", logicBorPoly_groove, logicMotherArc, kPhi, 10, th_tile, -5*th_tile );
	
    //To keep track of position along beam
    G4double layers_Xpos = 0.*cm;
    
    //Rotation matrix for Calorimeter
    G4RotationMatrix* r_y = new G4RotationMatrix();
    r_y->rotateY(90*deg);
    
    //Tiles angle of rotation w.r.t central Phi position
    G4double rotation_angle;
    
    //Tile index for copy number
    G4int tile_index = 0;
    
    //Indices to place correctly t0 where needed
    G4int t0_index = 0;
    G4double shift; //Y shift in planes without t0 w.r.t. rotated mother frame
   
    // edited by F. B.
    
	//Index to run on vector with tiles of random thickness but same type (see column definition of the matrix below)
	//For quadruplet
    int index = 0;
	//For triplet
	int index_trip = 0;
    
    //Matrices with t0->t3 LVolumes and t4->t6 LVolumes respectively (row: tile type; column: every entry is a tile of the same type but random thickness)
    std::vector<std::vector<G4LogicalVolume*>> logic_vec_t0_3;
    logic_vec_t0_3.push_back(logict0);
    logic_vec_t0_3.push_back(logict1);
    logic_vec_t0_3.push_back(logict2);
    logic_vec_t0_3.push_back(logict3);
    
    std::vector<std::vector<G4LogicalVolume*>> logic_vec_t4_6;
    logic_vec_t4_6.push_back(logict4);
    logic_vec_t4_6.push_back(logict5);
    logic_vec_t4_6.push_back(logict6);
    
    //Vector for LV fibers
    std::vector<G4LogicalVolume*> lv_fib_quad;
    lv_fib_quad.push_back(logicWLSGap0);
    lv_fib_quad.push_back(logicWLSGap1);
    lv_fib_quad.push_back(logicWLSGap2);
    lv_fib_quad.push_back(logicWLSGap3);
    
    std::vector<G4LogicalVolume*> lv_fib_trip;
    lv_fib_trip.push_back(logicWLSGap1_trip);
    lv_fib_trip.push_back(logicWLSGap2_trip);
    lv_fib_trip.push_back(logicWLSGap3_trip);
    
    
    //*** Optional offset in (Y, Z) plane ***
    G4double offset = 0.*mm;
	
  nLayers = 75; // <-------------------------------------------------------------- NLAYERS SELECTION!!!!!
//	Phi_layers = 1;
    //**** PHYSICAL PLACEMENTS ****
    for (G4int i = 0; i<nLayers; i++){
        if (i!=0) layers_Xpos += (layers_dist+zAbs);
        physiIronAbs = new G4PVPlacement(r_y, G4ThreeVector(layers_Xpos+zAbs/2.,demo_height,0.),logicIronAbs,"PhysiIronAbs",logicWorld,false,
										 2+i, fCheckOverlaps); //copy no range: [2, 76]
                                    
        physiBorPoly = new G4PVPlacement(r_y,G4ThreeVector(layers_Xpos+zAbs/2.,demo_height,0.),logicBorPoly,"PhysiBorPoly",logicWorld,false,
										 77+i, fCheckOverlaps); //copy no range: [77, 151]
		
		physiArcL = new G4PVPlacement(r_y, G4ThreeVector(zAbs+layers_dist/2.+layers_Xpos,demo_height,0.), logicArcL,"PhysiArcL",logicWorld, false,
									  0, fCheckOverlaps); //COPY NUMBER TO IMPROVE
		
		physiArcR = new G4PVPlacement(r_y, G4ThreeVector(zAbs+layers_dist/2.+layers_Xpos,demo_height,0.), logicArcR,"PhysiArcR",logicWorld, false,
									  0, fCheckOverlaps); //COPY NUMBER TO IMPROVE
         
        if (t0_index == 0 || t0_index == 1) {
            R_layers = 4;
            shift = 0*cm;
        } else if (t0_index>=2 && t0_index <= 4) {
            R_layers = 3;
            shift = h;
         }
		
        // ******* Loops to cycle on Phi-axis and R-axis to place physical volumes of tiles and fibers *******
        
        for (G4int j=0; j < Phi_layers; j++) {
			
			rotation_angle = th_tile*(4.5-j);
			G4RotationMatrix* r_x = new G4RotationMatrix();
			r_x->rotateX(-rotation_angle+pi/2.);
		
            for (G4int k=0; k < R_layers; k++) {
                
				if (i >= 5*N_z_instr) break; //Demo instrumentated only until 40th layer
				
                //To rotate Y-Z coordinates of tiles and fibers (called respectively x and y in TVector2)
                //Quadruplet and triplet fibers
                TVector2 old_ref, left_fib_t0, right_fib_t0, left_fib_t1, right_fib_t1, left_fib_t2, right_fib_t2, left_fib_t3, right_fib_t3;
                old_ref.Set(R_t0+h*(0.5+k)+shift+offset, 0.*cm+offset);
                left_fib_t0.Set( old_ref.X()+0.5*(fiberLenght0-h) , old_ref.Y()-0.5*fiberDistance );
                right_fib_t0.Set( old_ref.X()+0.5*(fiberLenght0-h) , old_ref.Y()+0.5*fiberDistance );
                left_fib_t1.Set( old_ref.X()+0.5*(fiberLenght1-h) , old_ref.Y()-0.5*fiberDistance );
                right_fib_t1.Set( old_ref.X()+0.5*(fiberLenght1-h) , old_ref.Y()+0.5*fiberDistance );
                left_fib_t2.Set( old_ref.X()+0.5*(fiberLenght2-h) , old_ref.Y()-0.5*fiberDistance+2*mm );
                right_fib_t2.Set( old_ref.X()+0.5*(fiberLenght2-h) , old_ref.Y()+0.5*fiberDistance+2*mm );
                left_fib_t3.Set( old_ref.X()+0.5*(fiberLenght3-h) , old_ref.Y()-0.5*fiberDistance-2*mm );
                right_fib_t3.Set( old_ref.X()+0.5*(fiberLenght3-h) , old_ref.Y()+0.5*fiberDistance-2*mm );
                //Implementing Vectors to save t4,t5,t6 fibers coordinates
                TVector2 left_fib_t4, right_fib_t4, left_fib_t5, right_fib_t5, left_fib_t6, right_fib_t6;
                left_fib_t4.Set( old_ref.X()+0.5*(fiberLenght1-h) , old_ref.Y()-0.5*fiberDistance );
                right_fib_t4.Set( old_ref.X()+0.5*(fiberLenght1-h) , old_ref.Y()+0.5*fiberDistance );
                left_fib_t5.Set( old_ref.X()+0.5*(fiberLenght2-h) , old_ref.Y()-0.5*fiberDistance );
                right_fib_t5.Set( old_ref.X()+0.5*(fiberLenght2-h) , old_ref.Y()+0.5*fiberDistance );
                left_fib_t6.Set( old_ref.X()+0.5*(fiberLenght3-h) , old_ref.Y()-0.5*fiberDistance+2*mm );
                right_fib_t6.Set( old_ref.X()+0.5*(fiberLenght3-h) , old_ref.Y()+0.5*fiberDistance+2*mm );
                
                TVector2 new_ref, Rleft_fib_t0, Rright_fib_t0, Rleft_fib_t1, Rright_fib_t1, Rleft_fib_t2, Rright_fib_t2, Rleft_fib_t3, Rright_fib_t3;
                TVector2 Rleft_fib_t4, Rright_fib_t4, Rleft_fib_t5, Rright_fib_t5, Rleft_fib_t6, Rright_fib_t6;
                new_ref = old_ref.Rotate(rotation_angle);
                Rleft_fib_t0 = left_fib_t0.Rotate(rotation_angle);
                Rright_fib_t0 = right_fib_t0.Rotate(rotation_angle);
                Rleft_fib_t1 = left_fib_t1.Rotate(rotation_angle);
                Rright_fib_t1 = right_fib_t1.Rotate(rotation_angle);
                Rleft_fib_t2 = left_fib_t2.Rotate(rotation_angle);
                Rright_fib_t2 = right_fib_t2.Rotate(rotation_angle);
                Rleft_fib_t3 = left_fib_t3.Rotate(rotation_angle);
                Rright_fib_t3 = right_fib_t3.Rotate(rotation_angle);
                
                Rleft_fib_t4 = left_fib_t4.Rotate(rotation_angle);
                Rright_fib_t4 = right_fib_t4.Rotate(rotation_angle);
                Rleft_fib_t5 = left_fib_t5.Rotate(rotation_angle);
                Rright_fib_t5 = right_fib_t5.Rotate(rotation_angle);
                Rleft_fib_t6 = left_fib_t6.Rotate(rotation_angle);
                Rright_fib_t6 = right_fib_t6.Rotate(rotation_angle);
				
                //***************** QUADRUPLET *****************
                if (R_layers==4) {
					//Now we need to find the maximum thickness value of the tiles of the same quadruplet, because X-coordinate could be higher or lower (due to randomization process of tiles thickness). Remember that tiles are placed in front of absorber layer, so, if a tile is thicker than the one below in the same quadruplet, the quadruplet physical placement needs to take into account maximum thickness value to place tiles along X-axis in order to avoid overlaps between tiles and absorber. This is also useful to place fibers always in the same position (along X) in the same quadruplet to go straight through the grooves of tiles above them.
					G4double max_tile_thickness=0.*mm;
					for(int it = 0; it<R_layers; it++){ //To find maximum thickness value of 4 type of tiles in the same quadruplet ("index" here is fixed)
						if (max_tile_thickness< random_thickness[it][index]) {
							max_tile_thickness = random_thickness[it][index];
						}
					}
					
                    physiTile = new G4PVPlacement(r_x,
                                                  G4ThreeVector(layers_Xpos+zAbs+ max_tile_thickness/2 , new_ref.X()+demo_height, new_ref.Y() ),
                                                  logic_vec_t0_3[k][index],
                                                  "PhysiTile",
                                                  logicWorld,
                                                  false,
                                                  1000+tile_index,  //copy no
												  fCheckOverlaps);
                    
                    // link copy no with tile coordinates: useful when collecting energy deposits.
                    // For quadruplets: assign R = -1 to t0 and (0, 1, 2) to LCMs
                    std::vector<G4int> v_coord{k-1, j, i};
                    gapCoord[1000 + tile_index] = v_coord;

                    //-------------------------------- FIBER PHYSICAL PLACEMENTS IN QUADRUPLET ---------------------------
                    // *** T0 FIBERS ***
                    if (k == 0) {
                        physiWLSfib = new G4PVPlacement(r_x, G4ThreeVector(layers_Xpos+zAbs+ max_tile_thickness -(grooveDiameter + abs(max_tile_thickness-tile_thickness_nom) )/2. , Rleft_fib_t0.X()+demo_height, Rleft_fib_t0.Y() ), lv_fib_quad[k],"PhysiFib_t0L",logicWorld,false, 0, fCheckOverlaps); //"0" copy number TO IMPROVE!!
                        
                        physiWLSfib = new G4PVPlacement(r_x, G4ThreeVector(layers_Xpos+zAbs+ max_tile_thickness -(grooveDiameter + abs(max_tile_thickness-tile_thickness_nom) )/2. , Rright_fib_t0.X()+demo_height, Rright_fib_t0.Y() ), lv_fib_quad[k],"PhysiFib_t0R",logicWorld,false, 0, fCheckOverlaps); //"0" copy number TO IMPROVE!!
                        
                    }

                    // *** T1 FIBERS ***
                    if (k == 1) {
                        physiWLSfib = new G4PVPlacement(r_x, G4ThreeVector(layers_Xpos+zAbs +(grooveDiameter + abs(max_tile_thickness-tile_thickness_nom) )/2. , Rleft_fib_t1.X()+demo_height, Rleft_fib_t1.Y() ), lv_fib_quad[k],"PhysiFib_t1L",logicWorld,false, 0, fCheckOverlaps); //"0" copy number TO IMPROVE!!
                        
                        physiWLSfib = new G4PVPlacement(r_x, G4ThreeVector(layers_Xpos+zAbs +(grooveDiameter + abs(max_tile_thickness-tile_thickness_nom) )/2. , Rright_fib_t1.X()+demo_height , Rright_fib_t1.Y() ), lv_fib_quad[k],"PhysiFib_t1R",logicWorld,false, 0, fCheckOverlaps); //"0" copy number TO IMPROVE!!
                        
                    }
                    // *** T2 FIBERS ***
                    if (k == 2) {
                        physiWLSfib = new G4PVPlacement(r_x, G4ThreeVector(layers_Xpos+zAbs+ max_tile_thickness -(grooveDiameter + abs(max_tile_thickness-tile_thickness_nom) )/2. , Rleft_fib_t2.X()+demo_height, Rleft_fib_t2.Y() ), lv_fib_quad[k],"PhysiFib_t2L",logicWorld,false, 0, fCheckOverlaps); //"0" copy number TO IMPROVE!!
                        
                        physiWLSfib = new G4PVPlacement(r_x, G4ThreeVector(layers_Xpos+zAbs+ max_tile_thickness -(grooveDiameter + abs(max_tile_thickness-tile_thickness_nom) )/2. , Rright_fib_t2.X()+demo_height, Rright_fib_t2.Y() ), lv_fib_quad[k],"PhysiFib_t2R",logicWorld,false, 0, fCheckOverlaps); //"0" copy number TO IMPROVE!!
                        
                    }
                    // *** T3 FIBERS ***
                    if (k == 3) {
                        physiWLSfib = new G4PVPlacement(r_x, G4ThreeVector(layers_Xpos+zAbs +(grooveDiameter + abs(max_tile_thickness-tile_thickness_nom) )/2. , Rleft_fib_t3.X()+demo_height, Rleft_fib_t3.Y() ), lv_fib_quad[k],"PhysiFib_t3L",logicWorld,false, 0, fCheckOverlaps); //"0" copy number TO IMPROVE!!
                        
                        physiWLSfib = new G4PVPlacement(r_x, G4ThreeVector(layers_Xpos+zAbs +(grooveDiameter + abs(max_tile_thickness-tile_thickness_nom) )/2. , Rright_fib_t3.X()+demo_height, Rright_fib_t3.Y() ), lv_fib_quad[k],"PhysiFib_t3R",logicWorld,false, 0, fCheckOverlaps); //"0" copy number TO IMPROVE!!
						
//						G4cout << "Max tile thick : "<< max_tile_thickness/mm << G4endl; //DEBUG
						
						//**** Creating a single module of BPE with single fiber holes for quadruplet layers ****
						
//						DetectorConstruction::BPEModule_holes(max_tile_thickness/*+1*um*/); //"BooleanProcessor::create Polyhedron:too many edges" error if you just give max_tile_thickness value. No visualization at all... 1*um doesn't affect positioning (just 3 orders of magnitude less)
//
//						//Logical and physical volumes of a single module of BPE with holes
//						logicBorPoly_groove = new G4LogicalVolume(solidBorPoly_groove, borpolyMaterial, "LogicBoratedPolyethylene");
//
//						logicBorPoly_groove->SetVisAttributes(borpolyColour);
//
//						G4RotationMatrix* r_module = new G4RotationMatrix();
//						r_module->rotateX(-rotation_angle);
//						r_module->rotateY(-90*deg);
						// Old BPE module physical placement
//						physiBorPoly_groove = new G4PVPlacement(r_module,G4ThreeVector(zAbs+layers_dist/2.+layers_Xpos,demo_height,0.), logicBorPoly_groove,
//																"PhysiBorPoly_module", logicWorld,false,
//																152,
//																fCheckOverlaps); //copy no range: [152, ...] TO IMPROVE!!!!
						
						
                        
                    }
                    
                }
                //***************** TRIPLET ******************
                else if (R_layers==3) {
					//Again we are looking for maximum tile thickness in the triplet (see comment for "QUADRUPLET" case)
					G4double max_tile_thickness=0.*mm;
					for(int it = 0; it<R_layers; it++){ //To find maximum thickness value of the tiles in the same quadruplet ("index_trip" here is fixed)
						if (max_tile_thickness< random_thickness_trip[it][index_trip]) {
							max_tile_thickness = random_thickness_trip[it][index_trip];
						}
					}
                    physiTile = new G4PVPlacement(r_x,
                                                  G4ThreeVector(layers_Xpos+zAbs+ max_tile_thickness/2. , new_ref.X()+demo_height, new_ref.Y() ),
                                                  logic_vec_t4_6[k][index_trip],
                                                  "PhysiTile",
                                                  logicWorld,
                                                  false,
                                                  1000+tile_index,  //copy no
												  fCheckOverlaps);

                    // link copy no with tile coordinates: useful when collecting energy deposits.
                    // For quadruplets: assign R = (0, 1, 2) to LCMs
                    std::vector<G4int> v_coord{k, j, i};
                    gapCoord[1000 + tile_index] = v_coord;

                    //-------------------------------- FIBER PHYSICAL PLACEMENTS IN TRIPLET ---------------------------
                    // *** T4 FIBERS ***
                    if (k == 0) {
                        physiWLSfib = new G4PVPlacement(r_x, G4ThreeVector(layers_Xpos+zAbs+ max_tile_thickness -(grooveDiameter + abs(max_tile_thickness-tile_thickness_nom) )/2. , Rleft_fib_t4.X()+demo_height, Rleft_fib_t4.Y() ), lv_fib_trip[k],"PhysiFib_t4L",logicWorld,false, 0, fCheckOverlaps); //"0" copy number TO IMPROVE!!
                        
                        physiWLSfib = new G4PVPlacement(r_x, G4ThreeVector(layers_Xpos+zAbs+ max_tile_thickness -(grooveDiameter + abs(max_tile_thickness-tile_thickness_nom) )/2. , Rright_fib_t4.X()+demo_height, Rright_fib_t4.Y() ), lv_fib_trip[k],"PhysiFib_t4R",logicWorld,false, 0, fCheckOverlaps); //"0" copy number TO IMPROVE!!
                        
                    }
                    // *** T5 FIBERS ***
                    if (k == 1) {
                        physiWLSfib = new G4PVPlacement(r_x, G4ThreeVector(layers_Xpos+zAbs+ (grooveDiameter + abs(max_tile_thickness-tile_thickness_nom) )/2. , Rleft_fib_t5.X()+demo_height, Rleft_fib_t5.Y() ), lv_fib_trip[k],"PhysiFib_t5L",logicWorld,false, 0, fCheckOverlaps); //"0" copy number TO IMPROVE!!
                        
                        physiWLSfib = new G4PVPlacement(r_x, G4ThreeVector(layers_Xpos+zAbs+ (grooveDiameter + abs(max_tile_thickness-tile_thickness_nom) )/2. , Rright_fib_t5.X()+demo_height, Rright_fib_t5.Y() ), lv_fib_trip[k],"PhysiFib_t5R",logicWorld,false, 0, fCheckOverlaps); //"0" copy number TO IMPROVE!!
                        
                    }
                    // *** T6 FIBERS ***
                    if (k == 2) {
                        physiWLSfib = new G4PVPlacement(r_x, G4ThreeVector(layers_Xpos+zAbs+ max_tile_thickness -(grooveDiameter + abs(max_tile_thickness-tile_thickness_nom) )/2. , Rleft_fib_t6.X()+demo_height, Rleft_fib_t6.Y() ), lv_fib_trip[k],"PhysiFib_t6L",logicWorld,false, 0, fCheckOverlaps); //"0" copy number TO IMPROVE!!

                        physiWLSfib = new G4PVPlacement(r_x, G4ThreeVector(layers_Xpos+zAbs+ max_tile_thickness -(grooveDiameter + abs(max_tile_thickness-tile_thickness_nom) )/2. , Rright_fib_t6.X()+demo_height, Rright_fib_t6.Y() ), lv_fib_trip[k],"PhysiFib_t6R",logicWorld,false, 0, fCheckOverlaps); //"0" copy number TO IMPROVE!!
						
//						//**** Creating a single module of BPE with single fiber holes for triplet layers ****
//						DetectorConstruction::BPEModule_holes(max_tile_thickness/*+2*um*/);
//
//						//Logical and physical volumes of a single module of BPE with holes
//						logicBorPoly_groove = new G4LogicalVolume(solidBorPoly_groove, borpolyMaterial, "LogicBoratedPolyethylene");
//
//						logicBorPoly_groove->SetVisAttributes(borpolyColour);
//
//						G4RotationMatrix* r_module = new G4RotationMatrix();
//						r_module->rotateX(-rotation_angle);
//						r_module->rotateY(-90*deg);
//
//						physiBorPoly_groove = new G4PVPlacement(r_module,G4ThreeVector(zAbs+layers_dist/2.+layers_Xpos,demo_height,0.), logicBorPoly_groove,
//																"PhysiBorPoly_module", logicWorld,false,
//																152,
//																fCheckOverlaps); //copy no range: [152, ...] TO IMPROVE!!!!

                    }
                    
                }
                
                
                
                tile_index ++; //for copy number only
                old_ref.Clear();
                new_ref.Clear();
				
            } //Closing R_layers loop
			
			//**** Creating a single module of BPE with two big holes for every Phi layer ****
			DetectorConstruction::BPEModule();

			//Logical and physical volumes of a single module of BPE with holes
			logicBorPoly_groove = new G4LogicalVolume(solidBorPoly_groove, borpolyMaterial, "LogicBoratedPolyethylene");

			logicBorPoly_groove->SetVisAttributes(borpolyColour);

			G4RotationMatrix* r_module = new G4RotationMatrix();
			r_module->rotateX(-rotation_angle);
			r_module->rotateY(-90*deg);

			physiBorPoly_groove = new G4PVPlacement (r_module,G4ThreeVector(zAbs+layers_dist/2.+layers_Xpos,demo_height,0.),
													 logicBorPoly_groove,"PhysiBorPoly_module", logicWorld,
													 false, 152,fCheckOverlaps); //COPY NO TO IMPROVE!!!
			
			
			
			//to move through the vector of tyles with random thickness
			if (t0_index == 0 || t0_index == 1) {
				index++;
			} else if (t0_index>=2 && t0_index <= 4) {
				index_trip++;
			 }
			
//			if (i >= 5*N_z_instr) {
//				//**** Creating a single module of BPE with single fiber holes for layers which are not instrumentated  ****
//				DetectorConstruction::BPEModule_holes(tile_thickness_nom);
//
//				//Logical and physical volumes of a single module of BPE with holes
//				logicBorPoly_groove = new G4LogicalVolume(solidBorPoly_groove, borpolyMaterial, "LogicBoratedPolyethylene");
//
//				logicBorPoly_groove->SetVisAttributes(borpolyColour);
//
//				G4RotationMatrix* r_module = new G4RotationMatrix();
//				r_module->rotateX(-rotation_angle);
//				r_module->rotateY(-90*deg);
				// OLD PHYSICAL PLACEMENT
//				physiBorPoly_groove = new G4PVPlacement(r_module,G4ThreeVector(zAbs+layers_dist/2.+layers_Xpos,demo_height,0.), logicBorPoly_groove,
//														"PhysiBorPoly_module_no_detector", logicWorld,false,
//														999,
//														fCheckOverlaps); //copy no range: [999, ...] TO IMPROVE!!!!
//			}
			
			
        } //Closing Phi_layers loop
        
        
        
        /*
        
        
        
        ////
        
       G4RotationMatrix* rm = new G4RotationMatrix(0*deg, 90*deg, 90*deg);
      // G4RotationMatrix* rm = new G4RotationMatrix();
           // rm->rotateY(90*deg);
           // rm->rotateZ(180*deg);
            
            
        physiGap[0] = new G4PVPlacement(rm,
            G4ThreeVector(-10*cm, -10*cm, -10*cm),
            logict1,
            "physiGap",
            logicMother,
            false,
            999,
            fCheckOverlaps);
            
        */    
            
            
            /*
            
            
            G4RotationMatrix* rotM = physiGap[0]->GetRotation();
		rotM->rotateX(90*deg);*/
		
		//*** Physical Placement of central BPE arc with BPE_module replicas with fiber channels ***
//		physiMotherArc = new G4PVPlacement(r_y, G4ThreeVector(zAbs+layers_dist/2.+layers_Xpos,demo_height,0.), logicMotherArc,"PhysiMotherArc",logicWorld, false,
//										   i+1/*, fCheckOverlaps*/); //COPY NUMBER TO IMPROVE and OVERLAPS!!
		
        t0_index ++;
        if (t0_index == 5) t0_index=0;
        
    } //Closing Z_layers loop
}


void DetectorConstruction::ConstructMuCatcher() {

	G4VisAttributes* muoAttr= new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	G4VisAttributes* brickAttr= new G4VisAttributes(G4Colour(0.5,0.5,0.5));
	G4VisAttributes* boxAttr= new G4VisAttributes(G4Colour(0.1,0.1,0.8));
	muoAttr->SetVisibility(true);
	brickAttr->SetVisibility(true);
	boxAttr->SetVisibility(true);
	
	// Build Muon Catcher
	
	// The external box is an aluminum frame
	G4Box* solidMuBox1 = new G4Box("solidMuonBox1", MuonBoxX1/2., MuonBoxY1/2.,MuonBoxZ1/2.);
    	G4Box* solidMuBox2 = new G4Box("solidMuonBox2", MuonBoxX2/2.,MuonBoxY2/2.,MuonBoxZ2/2.);
    	
    	G4Box* hollowMuBox1 = new G4Box("hollowMuBox1", (MuonBoxX1 - 2.*3.1*cm)/2., (MuonBoxY1 - 2*3.0*cm)/2., MuonActiveVolZ1/2.);
    	G4Box* hollowMuBox2 = new G4Box("hollowMuBox1", (MuonBoxX2 - 2.*1.7*cm)/2., MuonActiveVolY2/2., MuonActiveVolZ2/2.);
    	
  	solidMuonBox1 = new G4SubtractionSolid("Muon Catcher Box 1", solidMuBox1, hollowMuBox1, 0, G4ThreeVector(0,0,0));
    	solidMuonBox2 = new G4SubtractionSolid("Muon Catcher Box 2", solidMuBox2, hollowMuBox2, 0, G4ThreeVector(0,0,0));
    	

	// muoncatcher (active volume)

    	solidMuon1 = new G4Box("solidMuon1", MuonActiveVolX1/2., MuonActiveVolY1/2., MuonActiveVolZ1/2.);
    	solidMuon2 = new G4Box("solidMuon2", MuonActiveVolX2/2., MuonActiveVolY2/2., MuonActiveVolZ2/2.);
    	
    	// concrete brick
    	solidBrick = new G4Box("solidBrick", BrickSizeX/2., BrickSizeY/2., BrickSizeZ/2.);
    
    	// muoncatcher
    	logicMuonBox1 = new G4LogicalVolume(solidMuonBox1, boxMaterial, "logicMuonBox1"); // aluminium
    	logicMuonBox2 = new G4LogicalVolume(solidMuonBox2, boxMaterial, "logicMuonBox2"); // aluminium
    	
  	logicMuon1 = new G4LogicalVolume(solidMuon1, MCMaterial, "logicMuon1");
  	logicMuon2 = new G4LogicalVolume(solidMuon2, MCMaterial, "logicMuon2");
  	
  	logicBrick = new G4LogicalVolume(solidBrick, brickMaterial, "logicBrick");
	
	logicMuonBox1->SetVisAttributes(boxAttr);
	logicMuonBox2->SetVisAttributes(boxAttr);
	logicMuon1->SetVisAttributes(muoAttr);
  	logicMuon2->SetVisAttributes(muoAttr);
  	logicBrick->SetVisAttributes(brickAttr);

    	//G4int nMuCatch = 0; 
    	
    	
    	G4double demo_length = zBorPoly * nLayers; // 2.25 cm * 75
    	
    	//Rotation matrix for Calorimeter
	G4RotationMatrix* r_y = new G4RotationMatrix();
    	r_y -> rotateY(90*deg);
    
    	// position mu catcher

    	// active volume mu catcher 1
    	physiMuon1 =
    	new G4PVPlacement(r_y,
    	G4ThreeVector(demo_length + 91.5*cm + MuonBoxZ1, 136.5*cm + MuonActiveVolY1/2. + 3.0*cm,0),
    	logicMuon1,
    	"physiMuon1",
    	logicWorld,
    	false,
    	1,
    	fCheckOverlaps);
      	
	// active volume mu catcher 2
    	physiMuon2 =
      	new G4PVPlacement(r_y,
      	G4ThreeVector(demo_length + 91.5*cm + MuonBoxZ1 + MuonBoxZ2 + 2*cm /* 2 cm is a horizontal offset*/, 136.5*cm + MuonActiveVolY2/2. + (MuonBoxY2 - MuonActiveVolY2)/2.,  (MuonActiveVolX1 - MuonActiveVolX2)/2. ),
      	logicMuon2,
      	"physiMuon2",
      	logicWorld,
      	false,
      	2,
      	fCheckOverlaps);
      	
      	// external box mu catcher 1
    	physiMuonBox1 =
    	new G4PVPlacement(r_y,
    	G4ThreeVector(demo_length + 91.5*cm + MuonBoxZ1, 136.5*cm + MuonBoxY1/2., - (MuonBoxX1 - MuonActiveVolX1 - 2 * 3.1*cm)/2.),
    	logicMuonBox1,
    	"physiMuonBox1",
    	logicWorld,
    	false,
    	3,
    	fCheckOverlaps);
    	
    	// external box mu catcher 2
	physiMuonBox2 =
      	new G4PVPlacement(r_y,
      	G4ThreeVector(demo_length + 91.5*cm + MuonBoxZ1 + MuonBoxZ2 + 2*cm /* 2 cm is a horizontal offset*/, 136.5*cm + MuonBoxY2/2., - (MuonBoxX1 - MuonActiveVolX1 - 2 * MuonActiveVolX2 - 2.2*cm /*I don't know why 2.2cm, but i t works for alignment */)/2. ),
      	logicMuonBox2,
      	"physiMuonBox2",
      	logicWorld,
      	false,
      	4,
      	fCheckOverlaps);
    	
      	physiBrick =
  	new G4PVPlacement(r_y,
  	G4ThreeVector(demo_length + 91.5*cm + MuonBoxZ1 - BrickSizeZ, 136.5*cm + BrickSizeY/2.,  - (BrickSizeX - MuonActiveVolX1 - 2 * 3.1*cm)/2. ),
  	logicBrick,
  	"physiBrick",
  	logicWorld,
  	false,
  	5,
  	fCheckOverlaps);
  	
  	G4PVPlacement* physiBrick2 =
  	new G4PVPlacement(r_y,
  	G4ThreeVector(demo_length + 91.5*cm + MuonBoxZ1 + MuonBoxZ1 + MuonBoxZ2 + 1*cm /* 1 cm is a horizontal offset*/ + BrickSizeZ/2., 136.5*cm + BrickSizeY/2., - (BrickSizeX - MuonActiveVolX1 - 2 * 3.1*cm)/2.),
  	logicBrick,
  	"physiBrick",
  	logicWorld,
  	false,
  	6,
  	fCheckOverlaps);
}
