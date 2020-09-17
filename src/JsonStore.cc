#include "JsonStore.hh"

#include "G4ThreeVector.hh"

#include <iostream>
#include <fstream>
#include <exception>
#include "G4SystemOfUnits.hh"

#include "json.hpp"


using json = nlohmann::json;

JsonStore::JsonStore(const char* fileName)
{
	fDetMaterial = "Enriched";
	fFileName = "Src1_Enr";
	fDetGeo = "flatBEGe";

	fID = 0;
	fSrc = 1;

	fDeadLayerThickness = 0.02 * mm;
	fpLayerThickness = 0.001 * mm;
	fParticleE = 5.489 * MeV;

	fNPS = 100000;
}

JsonStore::~JsonStore() {}

bool JsonStore::Configure(const char* fileName)
{
	std::ifstream inputFile(fileName);
	json jsonfile;

	// check file exist 
	if (!inputFile.fail()) {
		inputFile >> jsonfile;
		inputFile.close();
	}
	else {
		return false;
	}

	// helper function to ensure key exist in json object ot prevent early exception
	auto keyExist = [](json j, std::string k) { return j.find(k) != j.end(); };

	for (auto& item : jsonfile.items())
	{
		try
		{
			if (item.key() == "FileName")
			{
				fFileName = item.value();
			}
			if (item.key() == "DetectorType")
			{
				fDetMaterial = item.value();
			}
			if (item.key() == "DetectorGeo")
			{
				fDetGeo = item.value();
			}
			else if (item.key() == "SourceID")
			{
				fSrc = item.value();
			}
			else if (item.key() == "SourceEnergy")
			{
				fParticleE = item.value();
			}
			else if (item.key() == "ID")
			{
				fID = item.value();
			}
			else if (item.key() == "DeadLayerThickness")
			{
				fDeadLayerThickness = item.value();
			}
			else if (item.key() == "pLayerThickness")
			{
				fpLayerThickness = item.value();
			}
			else if (item.key() == "Beam")
			{
				fNPS = item.value();
			}
			else // log wrong key 
			{
				G4cout << "WARNING: wrong key in configuration <" << item.key() << ">. " << G4endl;
			}
		}
		catch (json::out_of_range&) {
			G4cout << "WARNING:Some required keys are not set in: <" << item.key() << ">. " << G4endl;
		}
	}
	return true;
}

void JsonStore::PrintInfo() {
	G4cout << "######################### Simulation Info #########################" << G4endl;
	G4cout << "FileName=" << fFileName << G4endl;
	G4cout << "DetectorType=" << fDetMaterial << G4endl;
	G4cout << "DetectorGeo=" << fDetGeo << G4endl;
	G4cout << "SourceID=" << fSrc << G4endl;
	G4cout << "ParticleEnergy=" << fParticleE << G4endl;
	G4cout << "Beam=" << fNPS << G4endl;
	G4cout << "DeadLayerThickness=" << fDeadLayerThickness << G4endl;
	G4cout << "pLayerThickness=" << fpLayerThickness << G4endl;
	G4cout << "ID=" << fID << G4endl;
	G4cout << "###################################################################" << G4endl;
}
