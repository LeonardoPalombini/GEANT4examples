//source for particle generation: 511keV photons

#include "generator.hh"
#include "Randomize.hh"


MyPrimGen::MyPrimGen(double gammaE) {
	fParticleGun = new G4ParticleGun(1);

	G4ParticleTable* table = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = table->FindParticle("gamma");	
	G4ThreeVector pos(0., 0., 0.);
	//G4ThreeVector mom(0., 0., 1.);

	fParticleGun->SetParticlePosition(pos);
	//fParticleGun->SetParticleMomentumDirection(mom);
	fParticleGun->SetParticleEnergy(gammaE*MeV);
	//if(G4UniformRand() < 0.7) fParticleGun->SetParticleEnergy(4917 * keV);
	//else  fParticleGun->SetParticleEnergy(5117 * keV);
	fParticleGun->SetParticleDefinition(particle);
}
MyPrimGen::~MyPrimGen() {
	delete fParticleGun;
}

void MyPrimGen::GeneratePrimaries(G4Event* ev) {

	G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();

	//G4double th = acos(1. - 5.0391e-4 * G4UniformRand()); 3in
	//G4double th = acos(1. - 2.2471e-4 * G4UniformRand()); //2in
	//G4double phi = 2. * 3.1416 * G4UniformRand();
	G4ThreeVector mom(0., 0., 1.);
	//mom.setRThetaPhi(1., th, phi);

	fParticleGun->SetParticleMomentumDirection(mom);

	fParticleGun->GeneratePrimaryVertex(ev);
}