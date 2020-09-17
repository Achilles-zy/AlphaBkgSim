#include "B1StackingAction.hh"

#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4IonTable.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4HadronicProcessType.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1StackingAction::B1StackingAction()
	:G4UserStackingAction(),
	fFullChain(false),
	fCurrentParent(0),
	fCurrentSplitId(0),
	fRadioActiveSecID(0)
{
	// for HPGe detector, 10us time window
	fTimeWindow = 10 * us;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1StackingAction::~B1StackingAction()
{}


G4ClassificationOfNewTrack B1StackingAction::ClassifyNewTrack(const G4Track* track)
{
	G4double charge = track->GetDefinition()->GetPDGCharge();
	G4int ID = track->GetTrackID();
	G4int parentID = track->GetParentID();
	G4double trackTime = track->GetGlobalTime();

	// reset fCurrentParent to 0 for primary (ID=1)
	if (ID == 1) {
		fCurrentParent = 0;
		return fUrgent;
	}

	if (parentID == fRadioActiveSecID && trackTime > fTimeWindow&& parentID != 0) {
		if (fFullChain == false) {
			return fKill;
			//return fUrgent;
		}
		else {
			//G4cout << "wait" << G4endl;
			return fWaiting;
			//return fUrgent;
		}
	}

	// secondaries not from radioactivedecay
	G4int processType = track->GetCreatorProcess()->GetProcessSubType();
	//G4cout << track->GetCreatorProcess()->GetProcessName() << G4endl;
	if (processType != fRadioactiveDecay) {
		return fUrgent;
	}

	// reset global time after decay that exceeds detector time window
	if (trackTime > fTimeWindow) {
		// G4cout << "current track time (ns): " << trackTime << G4endl;
		const_cast<G4Track*>(track)->SetGlobalTime(0.0);
		// reload energy variables for new stages in the chain
		if (parentID != fCurrentParent) {
			fCurrentParent = parentID;
		}

	}
	// put radioactive secondaries into stacks
	//
	if (charge > 2.) {
		// if no fullchain kill track else put ions to waiting stack      
		fRadioActiveSecID = track->GetTrackID();
		return fUrgent;

	}
	else {
		// all other secondaries goes to normal stack
		return fUrgent;
	}
}