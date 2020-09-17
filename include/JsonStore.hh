/* JsonStore is modified to accommodate complex input file
	+ only moderate variables are analysed w/o check
	+ other variables are all optional
*/

#pragma once

#include <string>
#include <vector>
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

class JsonStore
{
public:
	JsonStore(const char* fileName);
	~JsonStore();


	std::string GetDetMaterial() { return fDetMaterial; }
	std::string GetFileName() { return fFileName; }
	std::string GetDetGeo() { return fDetGeo; }
	G4int GetID() { return fID; }
	G4int GetSrcType(){ return fSrc; }
	G4double GetParticleEnergy() { return fParticleE; }
	G4double GetDeadLayerThickness() { return fDeadLayerThickness; }
	G4double GetpLayerThickness() { return fpLayerThickness; }
	G4long GetNPS() { return fNPS; }
	G4bool Configure(const char*);

	void PrintInfo();
	// Particle distribution method

//private:
	std::string fDetMaterial;
	std::string fDetGeo;
	std::string fFileName;

	G4int fID;
	G4int fSrc;

	G4double fDeadLayerThickness;
	G4double fpLayerThickness;
	G4double fParticleE;

	long fNPS;

}; 
