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
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "B1EventAction.hh"
#include "B1RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "g4root.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction(B1RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.),
  bulkdepcount(0),
  fID(0),
  outputdata(0)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event*)
{    
  fEdep = 0.;
  bulkdepcount = 0;
  outputdata.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event* event)
{   
  // accumulate statistics in run action
  //fRunAction->AddEdep(fEdep);
	G4int eventid = event->GetEventID();
	if (eventid % 1000000 == 0) {
		G4cout << eventid << G4endl;
	}
  if (bulkdepcount == 0) fRunAction->CountDeadLayerEvent();
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillH1(0, fEdep);
  analysisManager->FillH1(1, fEdep);
  analysisManager->FillH1(2, fEdep);
  analysisManager->FillH1(3, fEdep);
  //G4cout << "E=" << fEdep << G4endl;
  if (fEdep > 5.304 * MeV) {
	  //G4cout << "E=" << fEdep << G4endl;
	  for (int i = 0; i < outputdata.size(); i++)
	  {
		  //G4cout << outputdata[i] << G4endl;
	  }
  }

  if (fEdep > 2.4 * MeV && fEdep < 2.6 * MeV) {
	  //G4cout << "E=" << fEdep << G4endl;
	  for (int i = 0; i < outputdata.size(); i++)
	  {
		  //G4cout << outputdata[i] << G4endl;
	  }
  }
  //analysisManager->FillNtupleIColumn(1, 5);
  //G4H1* h = G4AnalysisReader::Instance()->GetH1(fID);
  //h->Fill(fEdep);
  //G4cout << fEdep << G4endl;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
