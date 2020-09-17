
#ifndef B1StackingAction_h
#define B1StackingAction_h 1

#include "G4UserStackingAction.hh"
#include "globals.hh"
#include "G4Track.hh"

class G4Track;
class B1StackingAction : public G4UserStackingAction 
{
public:
	B1StackingAction();
	virtual ~B1StackingAction();

	virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);
	void SetTimeWindow(G4double tWin) { fTimeWindow = tWin; }

private:
	//G4bool   fBiasingFlag;
	G4bool   fFullChain;
	G4int    fCurrentParent;
	G4int    fRadioActiveSecID;
	G4int    fCurrentSplitId;
	G4double fTimeWindow;
};

#endif