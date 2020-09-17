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
/// \file exampleB1.cc
/// \brief Main program of the B1 example

#include "B1DetectorConstruction.hh"
#include "B1ActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "QBBC.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "B1PhysicsList.hh"

#include "Randomize.hh"
#include "JsonStore.hh"
#include "g4root.hh"
#include "FTFP_BERT_HP.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4GenericBiasingPhysics.hh"
#include "QBBC.hh"
#include "G4CoulombScattering.hh"
#include "G4IonElasticPhysics.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  const char* fileName = argv[1];
  auto input = new JsonStore(fileName);
  input->Configure(fileName);
  input->PrintInfo();

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  //
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // Set mandatory initialization classes
  //
  // Detector construction
  auto detectorconstruction = new B1DetectorConstruction();
  detectorconstruction->SetDeadlayerThickness(input->GetDeadLayerThickness());
  detectorconstruction->SetpLayerThickness(input->GetpLayerThickness());
  detectorconstruction->SetMatType(input->GetDetMaterial());
  detectorconstruction->SetDetGeo(input->GetDetGeo());
  detectorconstruction->Construct();
  runManager->SetUserInitialization(detectorconstruction);

  // Physics list
  G4int ModelID = 2;
  if (ModelID == 1) {
	  G4VModularPhysicsList* physicsList = new B1PhysicsList();
	  physicsList->SetVerboseLevel(0);
	  runManager->SetUserInitialization(physicsList);
  }

  if (ModelID == 2) {
	  G4VModularPhysicsList* physicslist1 = new QBBC;
	  physicslist1->SetDefaultCutValue(0 * um);
	  physicslist1->SetCutValue(50 * nm, "gamma");
	  physicslist1->SetCutValue(50 * nm, "e-");
	  physicslist1->SetCutValue(50 * nm, "e+");
	  physicslist1->RegisterPhysics(new G4RadioactiveDecayPhysics());
	  physicslist1->ReplacePhysics(new G4EmLivermorePhysics());
	  physicslist1->RegisterPhysics(new G4EmExtraPhysics());
	  physicslist1->RegisterPhysics(new G4IonElasticPhysics());
	  physicslist1->SetVerboseLevel(0);
	  runManager->SetUserInitialization(physicslist1);
  }
 
  if (ModelID == 3) {
	  G4VModularPhysicsList* physicslist2 = new FTFP_BERT_HP;
	  physicslist2->RegisterPhysics(new G4RadioactiveDecayPhysics());
	  physicslist2->ReplacePhysics(new G4EmStandardPhysics_option4());
	  physicslist2->SetVerboseLevel(0);
	  runManager->SetUserInitialization(physicslist2);
  }




  //


  B1ActionInitialization* actioninitialization = new B1ActionInitialization(input);
  G4cout << "fID=" << input->GetID() << G4endl;
  G4cout << "fFileName=" << input->GetFileName() << G4endl;
  G4cout << "fEnergy=" << input->GetParticleEnergy() << G4endl;
  actioninitialization->ActConfig(input->GetID(), input->GetSrcType(), input->GetFileName(), input->GetParticleEnergy());
  actioninitialization->printInfo();
  runManager->SetUserInitialization(new B1ActionInitialization(input));

  //runManager->Initialize();
  
  // Initialize visualization
  //

  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  

  // Get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/run/verbose 0");
  UI->ApplyCommand("/event/verbose 0");
  UI->ApplyCommand("/tracking/verbose 0");
  UI->ApplyCommand("/control/execute init_vis.mac");
  // Construct the default run manager
  // Process macro or start UI session
  //
  //if ( !ui ) { 
  //  // batch mode
  //  G4String command = "/control/execute ";
  //  G4String fileName = argv[1];
  //  UI->ApplyCommand(command+fileName);
  //}
  G4VisManager* visManager = nullptr;
  if (ui) {
    // interactive mode
    //UI->ApplyCommand("/control/execute init_vis.mac");
	  visManager = new G4VisExecutive;
	  visManager->Initialize();
    ui->SessionStart();
    delete ui;
  }
  else {
	  runManager->BeamOn(input->GetNPS());
  }
  //runManager->BeamOn(input->GetNPS());
  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !
  
  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
