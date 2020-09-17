#include "B1PhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmExtraPhysics.hh"
#include "G4GenericBiasingPhysics.hh"
#include "QBBC.hh"

#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4RadioactiveDecay.hh"
#include "G4CoulombScattering.hh"
#include "FTFP_BERT.hh"
#include "G4IonTable.hh"
#include "G4NuclearStopping.hh"
#include "G4ParticleTable.hh"
#include "G4IonElasticPhysics.hh"

// guide for different physics lists
// http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/html/TrackingAndPhysics/physicsProcess.html

B1PhysicsList::B1PhysicsList() : G4VModularPhysicsList()
{
	SetVerboseLevel(0);
	// Default physics
	SetDefaultCutValue(0 * um);
	SetParticleCuts(0 * mm, "proton");
	RegisterPhysics(new G4DecayPhysics());
	// EM physics (low energy option 3)
	G4VPhysicsConstructor* EmPhysics = new G4EmStandardPhysics_option3();
	G4VProcess* CS = new G4CoulombScattering();
	RegisterPhysics(new G4EmStandardPhysics_option4());
	// Synchroton Radiation & GN Physics
	RegisterPhysics(new G4EmExtraPhysics());
	// Hadronic physics
	//RegisterPhysics(new FTFP_BERT_());
	RegisterPhysics(new G4HadronPhysicsFTFP_BERT_HP());

	// Hadronic Elastic
	RegisterPhysics(new G4HadronElasticPhysics());
	// Stopping and ion physics
	RegisterPhysics(new G4StoppingPhysics());
	RegisterPhysics(new G4IonPhysics());
	// Radio active decay physics and add user defined data base
	RegisterPhysics(new G4RadioactiveDecayPhysics());
	RegisterPhysics(new G4IonElasticPhysics());
	//RegisterProcess(CS, G4IonTable::GetIonTable()->GetIon(86, 222, 0));






	// Add biasing to physical list
	//G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics();
	//biasingPhysics->NonPhysicsBias("gamma");
	//RegisterPhysics(biasingPhysics);
	std::cout << "INFO: Construct the Physical List register !" << std::endl;
}

B1PhysicsList::~B1PhysicsList()
{
	std::cout << "INFO: Deconstruct the Physical List register !" << std::endl;
}
