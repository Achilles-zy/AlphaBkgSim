//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4IonTable.hh"

#include "G4SingleParticleSource.hh"
#include "G4SPSPosDistribution.hh"
#include "G4GeneralParticleSource.hh"
#include "math.h"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
SrcID:
1: Vertical upper center
2: Vertical side center
3: Angular-even upper center
4: Angular-even side center
*/

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), 
  fEnvelopeBox(0),
  fSrcID(1),
  fParticleEnergy(5.489 *MeV)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  fParticleSrc = new G4GeneralParticleSource();
  // default particle kinematic
  fSrcID = 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fParticleSrc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName = "alpha");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(fParticleEnergy);

  fParticleSrc->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
  fParticleSrc->GetCurrentSource()->GetEneDist()->SetMonoEnergy(fParticleEnergy);
  fParticleSrc->GetCurrentSource()->SetParticleDefinition(particle);



  //G4LogicalVolume* Src = G4LogicalVolumeStore::GetInstance()->GetVolume("Src");
  //if (envLV) fEnvelopeBox = dynamic_cast<G4Box*>(envLV->GetSolid());
  double GeChamfer = 2. * mm;
  double outerGeRadius = 31.35 * mm;
  double innerGeRadius = 1.50 * mm;
  double GeHeight1 = 60.80 * mm;
  double SrcThickness = 0.01 * mm;
  double lSmallValue = 0.01 * mm;
  double orbradius = 1.89 * mm;
  double innergrooveradius = 7.5 * mm;
  //G4cout << "SrcID" << fSrcID << G4endl;
  
  if (fSrcID == 1) {
	  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
	  fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, GeHeight1 / 2 + GeChamfer));
	  fParticleGun->GeneratePrimaryVertex(anEvent);
  }
  
  if (fSrcID == 2) {
	  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-1, 0., 0.));
	  fParticleGun->SetParticlePosition(G4ThreeVector(outerGeRadius, 0, 0));
	  fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  if (fSrcID == 3) {
	  G4double theta = twopi * G4UniformRand();
	  G4double phi = pi * G4UniformRand();
	  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(cos(phi) * cos(theta), cos(phi) * sin(theta), sin(phi)));
	  fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, GeHeight1 / 2 + GeChamfer));
	  fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  if (fSrcID == 4) {
	  G4double theta = twopi * G4UniformRand();
	  G4double phi = pi * G4UniformRand();
	  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(cos(phi) * cos(theta), cos(phi) * sin(theta), sin(phi)));
	  fParticleGun->SetParticlePosition(G4ThreeVector(outerGeRadius, 0, 0));
	  fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  if (fSrcID == 5) {
	  G4double phi1 = twopi * G4UniformRand();
	  G4double theta1 = pi * G4UniformRand();
	  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(cos(phi1) * sin(theta1), sin(phi1) * sin(theta1), cos(theta1)));
	  //G4cout << G4ThreeVector(cos(phi1) * sin(theta1), sin(phi1) * sin(theta1), cos(theta1)) << G4endl;
	  G4double phi2 = twopi * G4UniformRand();
	  G4double theta2 = pi / 2 * G4UniformRand();
	  fParticleGun->SetParticlePosition(G4ThreeVector(orbradius * cos(phi2) * sin(theta2), orbradius * sin(phi2) * sin(theta2), orbradius * cos(theta2) - GeHeight1 / 2 + 1 * mm));
	  fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  if (fSrcID == 6) {
	  G4double phi1 = twopi * G4UniformRand();
	  G4double theta1 = pi * G4UniformRand();
	  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(cos(phi1) * sin(theta1), sin(phi1) * sin(theta1), cos(theta1)));
	  G4double phi2 = twopi * G4UniformRand();
	  G4double theta2 = pi / 2 * G4UniformRand();
	  G4double radius = orbradius * G4UniformRand();
	  fParticleGun->SetParticlePosition(G4ThreeVector(radius * cos(phi2) * sin(theta2), radius * sin(phi2) * sin(theta2), radius * cos(theta2) - GeHeight1 / 2 + 1 * mm));
	  fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  //Ra226 Sphere
  if (fSrcID == 7) {
	  G4int Z = 88, A = 226;
	  G4double ionCharge = 0. * eplus;
	  G4double excitEnergy = 0. * keV;

	  G4ParticleDefinition* ion
		  = G4IonTable::GetIonTable()->GetIon(Z, A, excitEnergy);
	  fParticleGun->SetParticleDefinition(ion);
	  fParticleGun->SetParticleCharge(ionCharge);
	  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));
	  //G4cout << G4ThreeVector(cos(phi1) * sin(theta1), sin(phi1) * sin(theta1), cos(theta1)) << G4endl;
	  G4double phi2 = twopi * G4UniformRand();
	  G4double theta2 = pi / 2 * G4UniformRand();
	  fParticleGun->SetParticlePosition(G4ThreeVector(orbradius * cos(phi2) * sin(theta2), orbradius * sin(phi2) * sin(theta2), orbradius * cos(theta2) - GeHeight1 / 2 + 1 * mm));
	  fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  //Po210 Sphere
  if (fSrcID == 8) {
	  G4int Z = 84, A = 210;
	  G4double ionCharge = 0. * eplus;
	  G4double excitEnergy = 0. * keV;

	  G4ParticleDefinition* ion
		  = G4IonTable::GetIonTable()->GetIon(Z, A, excitEnergy);
	  fParticleGun->SetParticleDefinition(ion);
	  fParticleGun->SetParticleCharge(ionCharge);
	  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 0));
	  //G4cout << G4ThreeVector(cos(phi1) * sin(theta1), sin(phi1) * sin(theta1), cos(theta1)) << G4endl;
	  G4double theta = twopi * G4UniformRand();
	  G4double radius = innergrooveradius * G4UniformRand();
	  //fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, -GeHeight1 / 2));
	  fParticleGun->SetParticlePosition(G4ThreeVector(radius * cos(theta), radius * sin(theta), -GeHeight1 / 2));
	  fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  //Po210 flatBEGe
  if (fSrcID == 9) {
	  G4int Z = 84, A = 210;
	  G4double ionCharge = 0. * eplus;
	  G4double excitEnergy = 0. * keV;

	  G4ParticleDefinition* ion
		  = G4IonTable::GetIonTable()->GetIon(Z, A, excitEnergy);
	  fParticleGun->SetParticleDefinition(ion);
	  fParticleGun->SetParticleCharge(ionCharge);
	  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 0));
	  //G4cout << G4ThreeVector(cos(phi1) * sin(theta1), sin(phi1) * sin(theta1), cos(theta1)) << G4endl;

	  G4double theta = twopi * G4UniformRand();
	  G4double radius = innergrooveradius * G4UniformRand();
	  //fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, -GeHeight1 / 2));
	  fParticleGun->SetParticlePosition(G4ThreeVector(radius* cos(theta), radius* sin(theta), -GeHeight1 / 2));
	  fParticleGun->GeneratePrimaryVertex(anEvent);
  }
  
  if (fSrcID == 10) {
	  G4int Z = 84, A = 210;
	  G4double ionCharge = 0. * eplus;
	  G4double excitEnergy = 0. * keV;

	  G4ParticleDefinition* ion
		  = G4IonTable::GetIonTable()->GetIon(Z, A, excitEnergy);
	  fParticleGun->SetParticleDefinition(ion);
	  fParticleGun->SetParticleCharge(ionCharge);
	  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 0));
	  //G4cout << G4ThreeVector(cos(phi1) * sin(theta1), sin(phi1) * sin(theta1), cos(theta1)) << G4endl;
	  G4double theta = twopi * G4UniformRand();
	  G4double radius = innergrooveradius * G4UniformRand();
	  //fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, -GeHeight1 / 2));
	  fParticleGun->SetParticlePosition(G4ThreeVector(radius * cos(theta), radius * sin(theta), -GeHeight1 / 2));
	  fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  if (fSrcID == 11) {
	  G4int Z = 88, A = 226;
	  G4double ionCharge = 0. * eplus;
	  G4double excitEnergy = 0. * keV;

	  G4ParticleDefinition* ion
		  = G4IonTable::GetIonTable()->GetIon(Z, A, excitEnergy);
	  fParticleGun->SetParticleDefinition(ion);
	  fParticleGun->SetParticleCharge(ionCharge);
	  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 0));
	  //G4cout << G4ThreeVector(cos(phi1) * sin(theta1), sin(phi1) * sin(theta1), cos(theta1)) << G4endl;
	  G4double theta = twopi * G4UniformRand();
	  G4double radius = innergrooveradius * G4UniformRand();
	  //fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, -GeHeight1 / 2));
	  fParticleGun->SetParticlePosition(G4ThreeVector(radius * cos(theta), radius * sin(theta), -GeHeight1 / 2));
	  fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  if (fSrcID == 12) {
	  G4int Z = 86, A = 222;
	  G4double ionCharge = 0. * eplus;
	  G4double excitEnergy = 0. * keV;

	  G4ParticleDefinition* ion
		  = G4IonTable::GetIonTable()->GetIon(Z, A, excitEnergy);
	  fParticleGun->SetParticleDefinition(ion);
	  fParticleGun->SetParticleCharge(ionCharge);
	  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 0));
	  //G4cout << G4ThreeVector(cos(phi1) * sin(theta1), sin(phi1) * sin(theta1), cos(theta1)) << G4endl;
	  G4double theta = twopi * G4UniformRand();
	  G4double radius = innergrooveradius * G4UniformRand();
	  //fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, -GeHeight1 / 2));
	  fParticleGun->SetParticlePosition(G4ThreeVector(radius * cos(theta), radius * sin(theta), -GeHeight1 / 2));
	  fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  if (fSrcID == 13) {
	  G4int Z = 84, A = 218;
	  G4double ionCharge = 0. * eplus;
	  G4double excitEnergy = 0. * keV;

	  G4ParticleDefinition* ion
		  = G4IonTable::GetIonTable()->GetIon(Z, A, excitEnergy);
	  fParticleGun->SetParticleDefinition(ion);
	  fParticleGun->SetParticleCharge(ionCharge);
	  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 0));
	  //G4cout << G4ThreeVector(cos(phi1) * sin(theta1), sin(phi1) * sin(theta1), cos(theta1)) << G4endl;
	  G4double theta = twopi * G4UniformRand();
	  G4double radius = innergrooveradius * G4UniformRand();
	  //fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, -GeHeight1 / 2));
	  fParticleGun->SetParticlePosition(G4ThreeVector(radius * cos(theta), radius * sin(theta), -GeHeight1 / 2));
	  fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  if (fSrcID == 14) {
	  G4int Z = 82, A = 214;
	  G4double ionCharge = 0. * eplus;
	  G4double excitEnergy = 0. * keV;

	  G4ParticleDefinition* ion
		  = G4IonTable::GetIonTable()->GetIon(Z, A, excitEnergy);
	  fParticleGun->SetParticleDefinition(ion);
	  fParticleGun->SetParticleCharge(ionCharge);
	  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 0));
	  //G4cout << G4ThreeVector(cos(phi1) * sin(theta1), sin(phi1) * sin(theta1), cos(theta1)) << G4endl;
	  G4double theta = twopi * G4UniformRand();
	  G4double radius = innergrooveradius * G4UniformRand();
	  //fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, -GeHeight1 / 2));
	  fParticleGun->SetParticlePosition(G4ThreeVector(radius * cos(theta), radius * sin(theta), -GeHeight1 / 2));
	  fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  //Po210 flatBEGe
  if (fSrcID == 15) {
	  G4double phi1 = twopi * G4UniformRand();
	  G4double theta1 = pi * G4UniformRand();
	  fParticleSrc->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
	  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(cos(phi1) * sin(theta1), sin(phi1) * sin(theta1), cos(theta1)));
	  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));
	  //G4cout << G4ThreeVector(cos(phi1) * sin(theta1), sin(phi1) * sin(theta1), cos(theta1)) << G4endl;

	  G4double theta = twopi * G4UniformRand();
	  G4double radius = innergrooveradius * G4UniformRand();
	  fParticleSrc->SetParticlePosition(G4ThreeVector(0, 0, -GeHeight1 / 2 + 0.000001));
	  fParticleSrc->GetCurrentSource()->GetPosDist()->SetPosDisType("Plane");
	  fParticleSrc->GetCurrentSource()->GetPosDist()->SetPosDisShape("Circle");
	  fParticleSrc->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0, 0, -GeHeight1 / 2));
	  fParticleSrc->GetCurrentSource()->GetPosDist()->SetRadius(radius);
	  //fParticleGun->SetParticlePosition(G4ThreeVector(radius * cos(theta), radius * sin(theta), -GeHeight1 / 2));
	  fParticleSrc->GeneratePrimaryVertex(anEvent);
  }

  //fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

