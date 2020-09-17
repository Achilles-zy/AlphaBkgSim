#ifndef B1TrackingAction_h
#define B1TrackingAction_h 1
#include "G4UserTrackingAction.hh"
#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "B1EventAction.hh"
#include "globals.hh"


class B1RunAction;
class B1EventAction;

class B1TrackingAction : public G4UserTrackingAction
{
public:
	B1TrackingAction(B1EventAction* eventaction);
	~B1TrackingAction();
	virtual void PreUserTrackingAction(const G4Track*);
	virtual void PostUserTrackingAction(const G4Track*);

private:
	B1EventAction* fEventAction;
}; 
#endif