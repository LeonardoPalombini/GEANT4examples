#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  virtual ~EventAction();

  virtual void  BeginOfEventAction(const G4Event* event);
  virtual void    EndOfEventAction(const G4Event* event);

private:
  G4int fPrimaryPDG;
  G4double fPrimaryPosX;
  G4double fPrimaryPosY;
  G4double fPrimaryPosZ;
  G4double fPrimaryEne;
  G4double fPrimaryDirX;
  G4double fPrimaryDirY;
  G4double fPrimaryDirZ;
  G4double fEnergyAbsTot;
  G4double fEnergyGapTot;
  G4double fEnergyMu1;
  G4double fEnergyMu2; 
  G4double fEneSpessore; 
  G4int nbOfGaps;
  G4int nbOfAbsor;
  G4int nbOfFibers;
  G4int nbOfCalo;
  G4int absNum;
  G4int gapNum;
  G4int sPlaneNum;
  G4int sipmNum;
};

#endif
