// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02PrimaryGeneratorAction.cc,v 1.2 1999/12/15 14:49:21 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// 
// $ Id: $
//
// From exampleEmPhys2/MyPrimaryGeneratorAction.cc,v 1.1 1998/02/06  maire 

#include "ExN02PrimaryGeneratorAction.hh"

#include "ExN02DetectorConstruction.hh"
#include "ExN02PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"

ExN02PrimaryGeneratorAction::ExN02PrimaryGeneratorAction(ExN02DetectorConstruction* myDC)
:myDetector(myDC),rndmFlag("off")
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  gunMessenger = new ExN02PrimaryGeneratorMessenger(this);

// default particle

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle 
    = particleTable->FindParticle(particleName="proton");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  particleGun->SetParticleEnergy(3.0*GeV);
  G4double position = -0.5*(myDetector->GetWorldFullLength());
  particleGun->SetParticlePosition(G4ThreeVector(position,0.*cm,0.*cm));

}

ExN02PrimaryGeneratorAction::~ExN02PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

void ExN02PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  // The Target is a box placed at (0,0,0)
  //
  G4double x0 = -0.5*(myDetector->GetTargetFullLength());
  G4double y0 = 0.*cm, z0 = 0.*cm;
  if (rndmFlag == "on")
     {y0 = (myDetector->GetTargetFullLength())*(G4UniformRand()-0.5);
      z0 = (myDetector->GetTargetFullLength())*(G4UniformRand()-0.5);
     } 
  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  particleGun->GeneratePrimaryVertex(anEvent);
}

