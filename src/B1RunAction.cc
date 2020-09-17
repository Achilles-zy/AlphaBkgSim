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
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
// #include "B1Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>
#include <iostream>
#include <string> 
#include <stdlib.h>
#include <iomanip>
#include <g4root.hh>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction()
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.),
  fDeadLayerEventNb(0),
  fFilename("Src1_Enr"),
  fID(0)
{ 
	G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
	accumulableManager->RegisterAccumulable(fEdep);
	accumulableManager->RegisterAccumulable(fEdep2);
	accumulableManager->RegisterAccumulable(fDeadLayerEventNb);

	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->SetVerboseLevel(1);
	analysisManager->CreateH1("h0", "EnergySpectrum", 500, 1 * keV, 8000 * keV);
	analysisManager->CreateH1("h1", "EnergySpectrum", 500, 3500 * keV, 5500 * keV);
	analysisManager->CreateH1("h2", "EnergySpectrum", 500, 1 * keV, 8000 * keV);
	analysisManager->CreateH1("h3", "EnergySpectrum", 500, 3500 * keV, 8000 * keV);
	//analysisManager->SetFirstNtupleColumnId(0);
	//analysisManager->SetFirstNtupleId(0);
	//G4cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$NtupleID" << G4endl;
	//G4cout << analysisManager->GetFirstNtupleColumnId() << G4endl;
	//G4cout << analysisManager->GetFirstNtupleColumnId() << G4endl;
	analysisManager->CreateNtuple("Data", "SurfaceEventNb");
	analysisManager->CreateNtupleDColumn("DeadLayerThickness");
	analysisManager->CreateNtupleIColumn("SurfaceEventNb");
	//analysisManager->CreateNtupleIColumn("Test1");
	analysisManager->FinishNtuple();
	//analysisManager->FillNtupleDColumn(analysisManager->GetFirstNtupleColumnId(), analysisManager->GetFirstNtupleColumnId(), 5);
	
////////////////////////////////////////////////////////////////////////////////////////////////
	

// Register accumulable to the accumulable manager

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{
	delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->SetFileName(fFilename);
	analysisManager->OpenFile();

	//G4cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ID" << G4endl;
	G4cout << fFilename << G4endl;
	G4cout << fID << G4endl;
	G4cout << std::to_string(fID) << G4endl;
	//analysisManager->SetHistoDirectoryName("./histo/");

	
	//analysisManager->CreateH1("h0", "EnergySpectrum", 500, 1 * keV, 6000 * keV);
	//analysisManager->CreateH1("h1", "EnergySpectrum", 500, 1 * keV, 6000 * keV);
	//analysisManager->CreateH1("h2", "EnergySpectrum", 500, 1 * keV, 6000 * keV);
	//analysisManager->CreateH1("h3", "EnergySpectrum", 500, 1 * keV, 6000 * keV);
	//analysisManager->CreateH1("h4", "EnergySpectrum", 500, 1 * keV, 6000 * keV);
	//analysisManager->CreateH1("h5", "EnergySpectrum", 500, 1 * keV, 6000 * keV);
	//analysisManager->CreateH1("h6", "EnergySpectrum", 500, 1 * keV, 6000 * keV);
	//analysisManager->CreateH1("h7", "EnergySpectrum", 500, 1 * keV, 6000 * keV);




	/*
	std::ifstream f(fFilename+".root");
	G4cout << fFilename + ".root" << G4endl;
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	G4AnalysisReader* analysisReader = G4AnalysisReader::Instance();
	if (f.good() == false) {
		G4cout << "##################################################false" << G4endl;

		analysisManager->SetFileName(fFilename);
		analysisManager->OpenFile();
		analysisManager->SetVerboseLevel(0);
		analysisManager->CreateH1("h0", "EnergySpectrum", 500, 1 * keV, 6000 * keV);
		analysisManager->CreateH1("h0", "EnergySpectrum", 500, 1 * keV, 6000 * keV);
		analysisManager->CreateH1("h1", "EnergySpectrum", 500, 1 * keV, 6000 * keV);
		analysisManager->CreateH1("h2", "EnergySpectrum", 500, 1 * keV, 6000 * keV);
		analysisManager->CreateH1("h3", "EnergySpectrum", 500, 1 * keV, 6000 * keV);
		analysisManager->CreateH1("h4", "EnergySpectrum", 500, 1 * keV, 6000 * keV);
		analysisManager->CreateH1("h5", "EnergySpectrum", 500, 1 * keV, 6000 * keV);
		analysisManager->CreateH1("h6", "EnergySpectrum", 500, 1 * keV, 6000 * keV);
		analysisManager->CreateH1("h7", "EnergySpectrum", 500, 1 * keV, 6000 * keV);


		analysisManager->CreateNtuple("SurfaceEvent", "Data");
		analysisManager->CreateNtupleIColumn("SurfaceEventNb");
		analysisManager->CreateNtupleIColumn("DeadLayerThickness");
		analysisManager->FinishNtuple();
	}
	else {
		//analysisManager->OpenFile(fFilename);
		//analysisManager->GetH1Id("h0");
		G4cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$true" << G4endl;

		//G4H1* h0 = analysisReader->GetH1(analysisReader->ReadH1("h0", fFilename));
		//G4cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$test" << G4endl;
		////G4cout << analysisReader->ReadH1("h0", fFilename,"C:/Users/hp/Desktop/B1/linuxbuild") << G4endl;
		////G4H1* h1 = analysisReader->GetH1(analysisReader->ReadH1("h1", fFilename, "C:/Users/hp/Desktop/B1/linuxbuild"));
		//G4cout << analysisReader->ReadH1("h1", fFilename, "C:/Users/hp/Desktop/B1/linuxbuild") << G4endl;
		//G4cout << analysisReader->ReadH1("h2", fFilename, "C:/Users/hp/Desktop/B1/linuxbuild") << G4endl;
		//analysisManager->GetH1Id("h1");
		//G4H1* h2 = analysisReader->GetH1(analysisReader->ReadH1("h2", fFilename, "C:/Users/hp/Desktop/B1/linuxbuild"));
		//G4H1* h3 = analysisReader->GetH1(analysisReader->ReadH1("h3", fFilename, "C:/Users/hp/Desktop/B1/linuxbuild"));
		//G4H1* h4 = analysisReader->GetH1(analysisReader->ReadH1("h4", fFilename, "C:/Users/hp/Desktop/B1/linuxbuild"));
		//G4cout << analysisManager->GetH1Id("h1") << G4endl;
		//G4cout << analysisManager->GetH1Id("h2") << G4endl;
		analysisReader->ReadH1("h1", fFilename, "C:/Users/hp/Desktop/B1/linuxbuild");
		analysisReader->ReadH1("h2", fFilename, "C:/Users/hp/Desktop/B1/linuxbuild");
		analysisReader->ReadH1("h3", fFilename, "C:/Users/hp/Desktop/B1/linuxbuild");
		analysisReader->ReadH1("h4", fFilename, "C:/Users/hp/Desktop/B1/linuxbuild");
		analysisReader->ReadH1("h5", fFilename, "C:/Users/hp/Desktop/B1/linuxbuild");

		//analysisManager->OpenFile(fFilename);
		//G4cout << analysisManager->GetNofH1s()<< G4endl;
		//analysisManager->GetFirstH1Id();
		//analysisManager->SetFirstH1Id(0);
		//analysisReader->SetFirstH1Id(0);
		
	}*/


	// inform the runManager to save random number seed
	G4RunManager::GetRunManager()->SetRandomNumberStore(false);

	// reset accumulables to their initial values
	G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
	accumulableManager->Reset();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* aRun)
{
	G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
	accumulableManager->Merge();

	G4int DeadLayerEventNb = fDeadLayerEventNb.GetValue();
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	//analysisManager->CreateNtuple("Data", "SurfaceEventNb");
	//analysisManager->CreateNtupleDColumn("DeadLayerThickness");
	analysisManager->FillNtupleIColumn(1, DeadLayerEventNb);
	analysisManager->FillNtupleDColumn(0, fDeadLayerThickness);
	analysisManager->AddNtupleRow();

	analysisManager->Write();
	analysisManager->CloseFile();
	//G4cout << "############################################################" << G4endl;
	//G4cout << "DeadLayerEventNb=" << DeadLayerEventNb << G4endl;
	//G4cout << "############################################################" << G4endl;
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}

void B1RunAction::CountDeadLayerEvent()
{
	fDeadLayerEventNb += 1;
	//G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
	//G4cout << "DeadLayerEventNb=" << fDeadLayerEventNb.GetValue() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

