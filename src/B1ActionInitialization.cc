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
/// \file B1ActionInitialization.cc
/// \brief Implementation of the B1ActionInitialization class

#include "B1ActionInitialization.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1RunAction.hh"
#include "B1EventAction.hh"
#include "B1SteppingAction.hh"
#include "G4SystemOfUnits.hh"
#include "B1TrackingAction.hh"
#include "B1StackingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1ActionInitialization::B1ActionInitialization(JsonStore* in)
 : G4VUserActionInitialization(),
	fID(0),
	fSrc(1),
	fEnergy(5.489 * MeV),
	fFileName("test")
{
	input = in;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1ActionInitialization::~B1ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1ActionInitialization::BuildForMaster() const
{
  B1RunAction* runAction = new B1RunAction;
  runAction->SetFileName(input->GetFileName());
  runAction->SetID(input->GetID());
  runAction->SetDeadLayerThickness(input->GetDeadLayerThickness());
  
  SetUserAction(runAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1ActionInitialization::Build() const
{
	B1PrimaryGeneratorAction* primaryGenerator = new B1PrimaryGeneratorAction();
	primaryGenerator->setEnergy(input->GetParticleEnergy());
	primaryGenerator->setSrcID(input->GetSrcType());
	SetUserAction(primaryGenerator);

  B1RunAction* runAction = new B1RunAction;
  //G4cout << "buildFilename" << input->GetFileName() << G4endl;
  runAction->SetFileName(input->GetFileName());
  runAction->SetID(input->GetID());
  runAction->SetDeadLayerThickness(input->GetDeadLayerThickness());
  SetUserAction(runAction);
  
  B1EventAction* eventAction = new B1EventAction(runAction);
  //G4cout << "buildID" << input->GetID() << G4endl;
  eventAction->setID(input->GetID());
  SetUserAction(eventAction);
  B1TrackingAction* trackingaction = new B1TrackingAction(eventAction);
  SetUserAction(trackingaction);
  B1StackingAction* stackingaction = new B1StackingAction;
  SetUserAction(stackingaction);
  B1SteppingAction* steppingaction = new B1SteppingAction(eventAction);
  steppingaction->SetpLayerThickness(input->GetpLayerThickness());
  SetUserAction(steppingaction);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
