#include "B1TrackingAction.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4AccumulableManager.hh"

#include "G4Electron.hh"
#include"G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4Gamma.hh"
#include "G4Alpha.hh"

#include "G4IonTable.hh"
#include "G4VProcess.hh"
#include "G4VAuxiliaryTrackInformation.hh"
#include "g4root.hh"
#include "G4HadronicProcessType.hh"
B1TrackingAction::B1TrackingAction(B1EventAction* eventaction):
	G4UserTrackingAction(),
	fEventAction(eventaction)
{
}


B1TrackingAction::~B1TrackingAction()
{
}

void B1TrackingAction::PreUserTrackingAction(const G4Track* track)
{
	//G4cout << "Globlal:" << track->GetGlobalTime() << G4endl;
	//G4cout << "Local:" << track->GetLocalTime() << G4endl;
	//G4cout << "PreTrackingAction" << G4endl;
	//G4cout << "Process = " << track->GetCreatorProcess()->GetProcessName() << G4endl;
	//const G4VProcess* process = track->GetCreatorProcess();
	
	//G4String name = process->GetProcessName().c_str();
    //G4cout << "Process = " << name << G4endl;
	//G4cout << "PreTrackingActionMap" << G4endl;
	//std::map< G4int, G4VAuxiliaryTrackInformation* >* g4map = track->GetAuxiliaryTrackInformationMap();
	//std::map< G4int, G4VAuxiliaryTrackInformation* >::iterator it = g4map->begin();
	//G4int mapsize = g4map->size();
	//for (int i = 0; i < mapsize; i++) {

		//G4VAuxiliaryTrackInformation* info = g4map->at(i);
		//info->Print();
	//}
	//G4cout << track->GetAuxiliaryTrackInformationMap() << G4endl;
	//G4cout << "PreTrackingAction" << G4endl;
	//G4cout << "Process = " << track->GetCreatorProcess()->GetProcessName() << G4endl;
	//G4cout << "time=" << track->GetGlobalTime() / 1 * s << G4endl;
	G4double charge = track->GetDefinition()->GetPDGCharge();
	G4int ID = track->GetTrackID();
	G4int parentID = track->GetParentID();
	G4double trackTime = track->GetGlobalTime();
	G4int processType = track->GetCreatorProcess()->GetProcessSubType();
	// reset fCurrentParent to 0 for primary (ID=1)
	if (ID != 1) {
		if (processType == fRadioactiveDecay && charge > 2) {
			//fEventAction->AddEdep(track->GetKineticEnergy());

		}
	}




	// secondaries not from radioactivedecay

}

void B1TrackingAction::PostUserTrackingAction(const G4Track* track)
{
	//G4cout << "PostTrackingAction" << G4endl;
	//G4cout << "Process = " << track->GetCreatorProcess()->GetProcessName() << G4endl;
	//track->GetCreatorProcess()->GetProcessName();


}
