#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "TVector2.h"
#include "G4Types.hh"
#include "Randomize.hh"
#include <vector>

class G4Box;
class G4Tubs;
class G4Trd;
class G4SubtractionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4GlobalMagFieldMessenger;
//class G4RandGauss;

//const G4int maxTiles = 2000;
//const G4int maxFibers = 2000;
const G4int maxTiles = 200000;
const G4int maxFibers = 9;
const G4int maxMuon = 2;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  DetectorConstruction();
  virtual ~DetectorConstruction();
  virtual void ConstructSDandField();

  void SetAbsorThickness(G4double);
  void SetGapThickness(G4double);

  void SetWorldMaterial (const G4String&);
  void SetCalorSizeYZ   (G4double);

  G4VPhysicalVolume* Construct();

  void PrintCalorParameters();

  G4double GetWorldSizeX()           {return worldSizeX;};
  G4double GetWorldSizeYZ()          {return worldSizeYZ;};

  G4double GetCalorThickness()       {return calorX;};
  G4double GetCalorSizeYZ()          {return calorSizeYZ;};

  G4int GetNPhiInstr() { return N_phi_instr; };
  G4int GetNZInstr() { return N_z_instr; };
  G4int       GetNbOfAbsor()             {return nbOfAbsor;};
  G4int       GetNbOfGaps()              {return nbOfGaps;};
  G4int       GetNbOfPoli()              {return nbOfPoli;};
  G4int       GetNbOfXCalo()             {return nbOfXCalo;};
  G4int       GetNbOfYCalo()             {return nbOfYCalo;};
  G4int       GetNbOfZCalo()             {return nbOfZCalo;};
  G4int       GetNbOfXCaloH()             {return nbOfXCaloH;};
  G4int       GetNbOfYCaloH()             {return nbOfYCaloH;};
  G4int       GetNbOfZCaloH()             {return nbOfZCaloH;};
  G4int       GetNbOfFibers()            {return nbOfFibers;};
  G4int       GetGapNum()                {return gapNum;};
  G4int       GetAbsNum()                {return absNum;};
  G4int       GetPoliNum()                {return poliNum;};
  G4int       GetFibNum()                {return fibNum;};
  G4int       GetSPlaneNum()             {return sPlaneNum;};
  G4int       GetSipmNum()               {return sipmNum;};
  G4Material* GetAbsorMaterial()         {return absorMaterial;};
  G4double    GetAbsorThickness()        {return absX;};
  G4Material* GetGapMaterial()           {return gapMaterial;};
  G4double    GetGapThickness()          {return gapZ;};
  G4Material* GetSipmMaterial()          {return sipmMaterial;};
  std::map<G4int, std::vector<G4int>> GetGapCoord() const { return gapCoord;}

  const G4VPhysicalVolume* GetphysiWorld()      {return physiWorld;};
  const G4Material*        GetWorldMaterial()   {return worldMaterial;};
  const G4VPhysicalVolume* GetAbsorber(G4int i)        {return physiAbs[i];};
  const G4VPhysicalVolume* GetGap(G4int i)             {return physiGap[i];};

private:
  static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;

  G4bool             fCheckOverlaps;
  G4Material*        worldMaterial;
  G4Material*        absorMaterial;
  G4Material*        gapMaterial;
  G4Material*        boxMaterial;
  G4Material*        boxMaterial1;
  G4Material*        fiberMaterial;
  G4Material*        poliMaterial;
  G4Material*        cladding1Material;
  G4Material*        cladding2Material;
  G4Material*        sipmMaterial;
  G4Material*        sensorPlaneMaterial;
  G4Material*        airCushionMaterial;
  G4Material*        brickMaterial;
  G4Material*        MCMaterial;

  G4double AirGapThickness;
  G4double BoxSizeY;
  G4double BoxSizeZ;
  G4double BoxSizeY1;
  G4double BoxSizeZ1;
  G4double BoxSizeX;
  G4double BoxSizeX1;
  G4double BoxSizeX2;

  // muon catcher

  // external box
  G4double MuonBoxX1;
  G4double MuonBoxY1;
  G4double MuonBoxZ1;
  G4double MuonBoxX2;
  G4double MuonBoxY2;
  G4double MuonBoxZ2;
  // active volume
  G4double MuonActiveVolX1;
  G4double MuonActiveVolY1;
  G4double MuonActiveVolZ1;
  G4double MuonActiveVolX2;
  G4double MuonActiveVolY2;
  G4double MuonActiveVolZ2;

  // brick

  G4double BrickSizeX;
  G4double BrickSizeY;
  G4double BrickSizeZ;

  G4double FeX;
  G4double FeY;
  G4double FeZ;

  G4double           absorThickness;
  G4double           gapThickness;

  G4double           motherX;
  G4double           motherY;
  G4double           motherZ;

  G4double           absX;
  G4double           absY;
  G4double           absZ;
  G4double           gapX;
  G4double           gapY;
  G4double           gapZ;

  G4double           poliX;
  G4double           poliY;
  G4double           poliZ;

  G4double           grooveX;
  G4double           grooveY;
  G4double           grooveZ;
  G4int              nbOfPlanes;
  G4int              nbOfFibers;
  G4int              nbOfPoli;
  G4double           layerThickness;
  G4double           calorSizeYZ;
  G4double           calorSizeXZ;

  G4double 	     fiberLenght0;
  G4double 	     fiberLenght1;
  G4double 	     fiberLenght2;
  G4double 	     fiberLenght3;

  G4double           calorX;
  G4double           calorZ;
  G4double           sensPlaneThickness;
  G4double           totalCalorThickness;
  G4double           fiberDiameter;
  G4double           grooveDiameter;
  G4double           fiberDistance;
  G4int              fibersPerRow;
  G4double           fiberStart;

  G4double           worldSizeYZ;
  G4double           worldSizeX;

  // Geometrical values of Demonstrator components

  G4double          R_t0;
  G4double          th_tile;
  G4double          D_Rg;
  G4int             N_phi_instr;
  G4int             N_z_instr;
  G4int             N_LCM;
  G4int             N_tile_LCM = 1200;
  G4int             N_tile_t0 = 160;
  G4double          h;
  G4double          tile_thickness_nom;
  G4double          tile_thickness_eff;
  // t0
  G4double          t0_min_base;
  G4double          t0_maj_base;
  // t1
  G4double          t1_min_base;
  G4double          t1_maj_base;
  // t2
  G4double          t2_min_base;
  G4double          t2_maj_base;
  // t3
  G4double          t3_min_base;
  G4double          t3_maj_base;
  // t4
  G4double          t4_min_base;
  G4double          t4_maj_base;
  // t5
  G4double          t5_min_base;
  G4double          t5_maj_base;
  // t6
  G4double          t6_min_base;
  G4double          t6_maj_base;
  // absorbers
  G4double          Rint_abs;
  G4double          arch_angle_abs;
  G4double          D_Ra;
  G4double          abs_thickness;
  G4double          abs_distance;
  G4int             Ntot_abs;

  // BPE
  G4double          Rint_BPE;
  G4double          arch_angle_BPE;
  G4double          D_Rb;
  G4double          BPE_thickness;
  G4int             Ntot_BPE;

  G4Box*             solidWorld;
  G4LogicalVolume*   logicWorld;
  G4VPhysicalVolume* physiWorld;


  G4LogicalVolume*   logicBoxEsterna;
  G4VPhysicalVolume* physiBoxEsterna[3];

  G4Box*             solidMother;
  G4LogicalVolume*   logicMother;
  G4VPhysicalVolume* physiMother;

  G4Box*             solidAbs;
  G4LogicalVolume*   logicAbs;
  G4VPhysicalVolume* physiAbs[maxTiles];
  
  G4Box*             solidPoli; 
  G4LogicalVolume*   logicPoli;
  G4VPhysicalVolume* physiPoli[maxTiles];

  G4Box*             solidGap;

  // vectors of solids of tiles
  std::vector<G4SubtractionSolid*> t0;
  std::vector<G4SubtractionSolid*> t1;
  std::vector<G4SubtractionSolid*> t2;
  std::vector<G4SubtractionSolid*> t3;
  std::vector<G4SubtractionSolid*> t4;
  std::vector<G4SubtractionSolid*> t5;
  std::vector<G4SubtractionSolid*> t6;

  // vectors of logical volumes of tiles

  std::vector<G4LogicalVolume*>   logict0;
  std::vector<G4LogicalVolume*>   logict1;
  std::vector<G4LogicalVolume*>   logict2;
  std::vector<G4LogicalVolume*>   logict3;
  std::vector<G4LogicalVolume*>   logict4;
  std::vector<G4LogicalVolume*>   logict5;
  std::vector<G4LogicalVolume*>   logict6;

  G4VPhysicalVolume* physiGap[maxTiles];

  G4Tubs*            solidFiber;
  G4LogicalVolume*   logicFiber;
  G4VPhysicalVolume* physiFiber;

  G4Tubs*            solidClad1;
  G4LogicalVolume*   logicClad1;
  G4VPhysicalVolume* physiClad1;

  G4Tubs*            solidClad2;
  G4LogicalVolume*   logicClad2;
  G4VPhysicalVolume* physiClad2;

  G4Tubs*            solidHole;
  G4LogicalVolume*   logicHole;

  G4LogicalVolume*   logicWLSAbs;
  G4VPhysicalVolume* physiWLSAbs;

  G4LogicalVolume*   logicWLSGap0;
  G4LogicalVolume*   logicWLSGap1;
  G4LogicalVolume*   logicWLSGap2;
  G4LogicalVolume*   logicWLSGap3;
  // triplet
  G4LogicalVolume*   logicWLSGap0_trip;
  G4LogicalVolume*   logicWLSGap1_trip;
  G4LogicalVolume*   logicWLSGap2_trip;
  G4LogicalVolume*   logicWLSGap3_trip;
  G4VPhysicalVolume* physiWLSGap[maxFibers];
  G4VPhysicalVolume* physiWLSGap_trip[maxFibers];

  G4Box*             solidSipm;
  G4LogicalVolume*   logicSipm;
  G4VPhysicalVolume* physiSipm;

  G4Box*             solidAirCushion;
  G4LogicalVolume*   logicAirCushion;
  G4VPhysicalVolume* physiAirCushion;

  G4Box*             solidSensorPlane;
  G4LogicalVolume*   logicSensorPlane;
  G4VPhysicalVolume* physiSensorPlane[maxTiles];
  G4VPhysicalVolume* physiSensorPlaneDummy;

  G4Box*             solidBrick;
  G4LogicalVolume*   logicBrick;
  G4VPhysicalVolume* physiBrick;

  G4Box*             solidBox;
  G4LogicalVolume*   logicBox;
  G4VPhysicalVolume* physiBox;

  G4SubtractionSolid* solidMuonBox1;
  G4LogicalVolume* logicMuonBox1;
  G4VPhysicalVolume* physiMuonBox1;

  G4SubtractionSolid* solidMuonBox2;
  G4LogicalVolume* logicMuonBox2;
  G4VPhysicalVolume* physiMuonBox2;

  G4Box*             solidMuon1;
  G4LogicalVolume*   logicMuon1;
  G4VPhysicalVolume* physiMuon1;

  G4Box*             solidMuon2;
  G4LogicalVolume*   logicMuon2;
  G4VPhysicalVolume* physiMuon2;
  G4LogicalVolume*   logicSpessore;
  G4VPhysicalVolume* physiSpessore;

  G4Box*             solidFe;
  G4LogicalVolume*   logicFe;
  G4VPhysicalVolume* physiFe;

  G4LogicalVolume*   logicBoxFinal;
  G4VPhysicalVolume* physiBoxFinal[3];

  G4Box*             solidBox1;
  G4Box*             solidBox2;
  G4LogicalVolume*   logicBox1;
  G4VPhysicalVolume* physiBox1[3];

  G4Box*             solidMCScin;
  G4LogicalVolume*   logicMCScin;
  G4VPhysicalVolume* physiMCScin;

  G4int nbOfXCalo;
  G4int nbOfYCalo;
  G4int nbOfZCalo;
  G4int nbOfXCaloH;
  G4int nbOfYCaloH;
  G4int nbOfZCaloH;

  G4int nbOfGaps;
  G4int nbOfBoxes;
  G4int nbOfAbsor;
  G4double modulesGap ;
  G4int gapNum;
  G4int absNum;
  G4int poliNum; 
  G4int fibNum;
  G4int sPlaneNum;
  G4int sipmNum;
    
  //Iron absorber
  G4Material* ironAbsMaterial;
    
  //G4int R_t0;
  //G4int h;
  G4double arc_angle;
  G4double start_phi;
  //G4int D_Ra;
  G4double zAbs;
  G4double layers_dist;
  G4int nLayers;
    
  G4Tubs*             solidIronAbs;
  G4LogicalVolume*    logicIronAbs;
  G4VPhysicalVolume*  physiIronAbs;
    
  //Borated Polyehtylene arcs
  G4Material* borpolyMaterial;
    
  // G4double D_Rg;
  // G4double D_Rb;
  G4double zBorPoly;
    
  G4Tubs*             solidBorPoly;
  G4SubtractionSolid* solidBorPoly_groove;
  G4Tubs*				solidArcL;
  G4Tubs*				solidArcR;
  G4Tubs*				solidMotherArc;
  G4LogicalVolume*    logicBorPoly;
  G4LogicalVolume*    logicBorPoly_groove;
  G4LogicalVolume*	logicArcL;
  G4LogicalVolume*	logicArcR;
  G4LogicalVolume*	logicMotherArc;
  G4VPhysicalVolume*  physiBorPoly;
  G4VPhysicalVolume*  physiBorPoly_groove;
  G4VPhysicalVolume* 	physiArcL;
  G4VPhysicalVolume* 	physiArcR;
  G4VPhysicalVolume*	physiMotherArc;
    
  G4VPhysicalVolume* physiTile;
  G4VPhysicalVolume* physiWLSfib;
    
  G4int R_layers;
  G4int Phi_layers;

  // tile thickness randomization
  G4double sigma; // sigma for randomization
  G4double THICK_OFFSET;
  CLHEP::RandGauss rnd = CLHEP::RandGauss( *(CLHEP::HepRandom::getTheEngine()) , 0.0 , 1.0 );
  G4double demo_height; //height between ground and centre of the inner arc of iron absorber (along y-axis)
	
  //Vectors to save random thickness of tiles, useful for physical placements
  std::vector<std::vector<G4double>> random_thickness, random_thickness_trip;
  std::vector<G4double> t0_thick, t1_thick, t2_thick, t3_thick, t4_thick, t5_thick, t6_thick;

  std::map<G4int, std::vector<G4int>> gapCoord;
  //void BPEModule_holes(G4double);
  void BPEModule();
  void ConstructIronAbs_BorPoly();
  void ConstructBorPoly();

  void DefineMaterials();
  void ComputeCalorParameters();
  void ComputeDemonstratorParameters();
  void ConstructCalorimeter();
  void ConstructTiles();
  G4double RandomThickness(G4double mean);
  void ConstructMuCatcher();

  G4LogicalVolume* ConstructWLSFiber(G4double);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
