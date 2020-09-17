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
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4HadronicProcessType.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0),
	fBulk(0),
	fEnv(0),
	fpLayer(0),
	fpLayerThickness(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) { 
    const B1DetectorConstruction* detectorConstruction
      = static_cast<const B1DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fBulk = detectorConstruction->GetBulk();   
	fEnv = detectorConstruction->GetEnv();
	fpLayer = detectorConstruction->GetpLayer();

  }
  
  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();

      
  // check if we are in scoring volume
  //if (volume != fScoringVolume) return;

  // collect energy deposited in this step
  //G4double edepStep1 = step->GetTotalEnergyDeposit();
  //G4cout << "Edepstep=" << edepStep1 << G4endl;

  //G4double edepStep = step->GetTotalEnergyDeposit();
  //fEventAction->AddEdep(edepStep);
  //if (edepStep != 0) fEventAction->depCount();
  //if (step->GetTrack()->GetTrackID() == 2) {
	 // G4cout
		//  << "Particle = " << step->GetTrack()->GetParticleDefinition()->GetParticleName() << " "
		//  << "DeltaE = " << step->GetTotalEnergyDeposit() << " "
		//  << "Process = " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << " "
		//  << "Momentum Dir = " << step->GetPostStepPoint()->GetMomentumDirection() << " "
		//  << "Position = " << step->GetPostStepPoint()->GetPosition() << " "
		//  << "Volume = " << step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName()
		//  << G4endl;

  //}
  if (fpLayerThickness > 0) {
	  if (volume == fBulk) {
		  //if (volume != fEnv) {
		  G4double edepStep = step->GetTotalEnergyDeposit();
		  fEventAction->AddEdep(edepStep);

		  if (edepStep != 0) fEventAction->depCount();
	  }
  }
  if (fpLayerThickness == 0) {
	  if (volume == fBulk || volume == fpLayer) {
		  //if (volume != fEnv) {
		  G4double edepStep = step->GetTotalEnergyDeposit();
		  fEventAction->AddEdep(edepStep);

		  if (edepStep != 0) fEventAction->depCount();
	  }
  }

  G4int processtype = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessSubType();
  G4int parentID = step->GetTrack()->GetParentID();
  if (parentID == 1 && processtype == fRadioactiveDecay) {
	  //step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
	  //step->GetTrack()->SetTrackStatus(fPostponeToNextEvent);
  }

  
  //G4String creatorName;
  //G4int trackID = step->GetTrack()->GetTrackID();
  //if (trackID == 1)  creatorName = "primary";
  //else creatorName = step->GetTrack()->GetCreatorProcess()->GetProcessName();
  ////G4cout << "prepointprocess = " << creatorName << G4endl;
  //G4String name = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  /*
  G4Track* track = step->GetTrack();
  G4int trackID = track->GetTrackID();
  G4int parentID = track->GetParentID();
  G4int stepID = track->GetCurrentStepNumber();
  G4double globalTime = track->GetGlobalTime();
  G4String particalName = track->GetDefinition()->GetParticleName();// equal to GetDynamicParticle()
  G4String volumeName = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
  G4String procName = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  // G4String particalName2 = track->GetDynamicParticle()->GetDefinition()->GetParticleName();
  G4ThreeVector dir = step->GetPostStepPoint()->GetMomentumDirection();
  G4ThreeVector xyzTrack = track->GetPosition();// equal to PostStepPoint()
  G4ThreeVector xyzPost = step->GetPostStepPoint()->GetPosition();// string type
  G4double edepStep = step->GetTotalEnergyDeposit();
  G4double StepLen = step->GetStepLength();
  //
  G4String output = std::to_string(trackID) + " ";
  output += std::to_string(parentID) + " ";
  output += std::to_string(stepID) + " ";
  output += std::to_string(edepStep) + " ";
  output += particalName + " ";
  output += volumeName + " ";
  output += procName + " ";
  output += "(" + std::to_string(dir[0]) + "," + std::to_string(dir[1]) + "," + std::to_string(dir[2]) + ")"+ " ";
  output += "(" + std::to_string(xyzPost[0]) + "," + std::to_string(xyzPost[1]) + "," + std::to_string(xyzPost[2]) + ")" + " ";
  output += std::to_string(StepLen);

  fEventAction->store(output);
  */


}

void B1SteppingAction::PreUserSteppingAction(const G4Step* step) {


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

