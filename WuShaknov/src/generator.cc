//source for particle generation: uncomment the desired physical process

#include "generator.hh"
#include "Randomize.hh"

//full 22Na decay
/*
MyPrimGen::MyPrimGen() {
	fParticleGun = new G4ParticleGun(1);

	G4ParticleTable* table = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = table->FindParticle("geantino");	//fictitious particle to hold place
	G4ThreeVector pos(0., 0., 0.);
	G4ThreeVector mom(0., 0., 0.);

	fParticleGun->SetParticlePosition(pos);
	fParticleGun->SetParticleMomentumDirection(mom);
	fParticleGun->SetParticleMomentum(0.);
	fParticleGun->SetParticleDefinition(particle);
}
MyPrimGen::~MyPrimGen() {
	delete fParticleGun;
}

void MyPrimGen::GeneratePrimaries(G4Event* ev) {

	G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();

	if (particle == G4Geantino::Geantino()) {
		G4int Z = 11;
		G4int A = 22;
		G4double charge = 0. * eplus;
		G4double energy = 0. * keV;
		G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z, A, energy);
		fParticleGun->SetParticleDefinition(ion);
		fParticleGun->SetParticleCharge(charge);
	}

	fParticleGun->SetParticlePosition(G4ThreeVector((G4UniformRand()*0.4-0.2)*cm, (G4UniformRand() * 0.4 - 0.2)*cm,0.));
	

	fParticleGun->GeneratePrimaryVertex(ev);
	count++;
	if(count%10000 == 0) G4cout << count << G4endl;
}
*/

//collinear polarized 511keV photons

MyPrimGen::MyPrimGen() {
	fParticleGun = new G4ParticleGun(1);

	G4ParticleTable* table = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = table->FindParticle("gamma");
	//G4ThreeVector pos(0., 0., 0.);
	//fParticleGun->SetParticleMomentum(100.*keV);
	fParticleGun->SetParticleDefinition(particle);

}
MyPrimGen::~MyPrimGen() {
	delete fParticleGun;
}

void MyPrimGen::GeneratePrimaries(G4Event* ev) {

	G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();

	fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
	
	/*if (G4UniformRand()<0.5) th = acos(1. - 2. * 0.003 * G4UniformRand());
	else th = acos(1. - 2.*(0.997 + 0.003 * G4UniformRand()));
	phi = 2. * 3.1416 * G4UniformRand();*/
	
	th = acos(1. - 2. * 0.0021 * G4UniformRand());
	phi = 2. * 3.1416 * G4UniformRand();

	G4ThreeVector mom1(0.,0.,1.);
	mom1.setRThetaPhi(1.,th,phi);
	
	fParticleGun->SetParticleEnergy(511. * keV);
	
	G4ThreeVector pol1 = mom1.orthogonal();
	pol1.rotate(mom1, 2. * 3.1416 * G4UniformRand());
	G4ThreeVector pol2 = pol1;
	pol2.rotate(mom1, 90 * deg);



	fParticleGun->SetParticleMomentumDirection(mom1);
	fParticleGun->SetParticlePolarization(pol1);
	fParticleGun->GeneratePrimaryVertex(ev);

	fParticleGun->SetParticleMomentumDirection(-mom1);
	fParticleGun->SetParticlePolarization(pol2);
	fParticleGun->GeneratePrimaryVertex(ev);

	//G4cout << "Polarizations: " << pol1 << "  " << pol2 << G4endl;
	

	/*th = acos(1. - 2. * G4UniformRand());
	phi = 2. * 3.1416 * G4UniformRand();
	G4ThreeVector momEx(sin(th) * cos(phi), sin(th) * sin(phi), cos(th));

	fParticleGun->SetParticleEnergy(1275. * keV);
	fParticleGun->SetParticleMomentumDirection(momEx);
	fParticleGun->GeneratePrimaryVertex(ev);*/

	if (count % 10000 == 0) G4cout << count << G4endl;
	count++;
	
}


//single 511keV photon
/*
MyPrimGen::MyPrimGen() {
	fParticleGun = new G4ParticleGun(1);

	G4ParticleTable* table = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = table->FindParticle("gamma");	
	G4ThreeVector pos(0., 46.6*cm, 12.5*cm);
	G4ThreeVector mom(0., -1., 0.);

	fParticleGun->SetParticlePosition(pos);
	fParticleGun->SetParticleMomentumDirection(mom);
	fParticleGun->SetParticleMomentum(1400.*keV);
	fParticleGun->SetParticleDefinition(particle);
}
MyPrimGen::~MyPrimGen() {
	delete fParticleGun;
}

void MyPrimGen::GeneratePrimaries(G4Event* ev) {

	G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
	th = acos(1. - 2. * 0.0016 * G4UniformRand());
	phi = 2. * 3.1416 * G4UniformRand();

	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(th) * sin(phi),-cos(th), sin(th)*cos(phi)));


	fParticleGun->GeneratePrimaryVertex(ev);
}
*/