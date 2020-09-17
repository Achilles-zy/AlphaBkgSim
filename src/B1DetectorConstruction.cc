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
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Torus.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "g4root.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
	: G4VUserDetectorConstruction(),
	fScoringVolume(0),
	fSrc(0),
	fBulk(0),
	fDeadLayer(0),
	fInnerDeadLayer(0),
	fEnv(0),
	fDetDeadlayer(1 * mm),
	fpLayerThickness(0.001 * mm),
	fMatType("Enriched"),
	fDetGeo("flatBEGe"),
	fpLayer(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  G4cout << "########################## Construction Info ##########################" << G4endl;
  G4cout << "DeadLayerThickness=" << fDetDeadlayer << G4endl;
  G4cout << "pLayerThickness=" << fpLayerThickness << G4endl;
  G4cout << "detGeo=" << fDetGeo << G4endl;
  G4cout << "#######################################################################" << G4endl;
  
  if (fpLayerThickness == 0) {
	  SetpLayerThickness(0.001 * nm);
  }
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  //vacuum
  G4Material* Air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* vacuum = new G4Material("Vacuum ", 0 * g / cm3, 1,
	  kStateGas, 50 * kelvin, 0 * atmosphere);
  vacuum->AddMaterial(Air, 1.);

  //Enriched Ge
  G4int A, Z, nAtoms, nComponents;
  G4double moleMass, density, fractionmass;
  auto Ge70 = new G4Isotope("Ge70", Z = 32, A = 70, moleMass = 69.9 * g / mole);
  auto Ge72 = new G4Isotope("Ge72", Z = 32, A = 72, moleMass = 71.92 * g / mole);
  auto Ge73 = new G4Isotope("Ge73", Z = 32, A = 73, moleMass = 72.92 * g / mole);
  auto Ge74 = new G4Isotope("Ge74", Z = 32, A = 74, moleMass = 73.92 * g / mole);
  auto Ge76 = new G4Isotope("Ge76", Z = 32, A = 76, moleMass = 75.92 * g / mole);
  auto GeEn = new G4Element("enrichedGermanium", "GeEn", 5);
  GeEn->AddIsotope(Ge70, 0.001);
  GeEn->AddIsotope(Ge72, 0.001);
  GeEn->AddIsotope(Ge73, 0.001);
  GeEn->AddIsotope(Ge74, 0.130);
  GeEn->AddIsotope(Ge76, 0.867);
  auto matGeEn = new G4Material("GeEn", 5.539 * g / cm3, 1);
  matGeEn->AddElement(GeEn, 1.0);

  auto matGeNa = new G4Material("GeNa", 32, 72.61 * g / mole, 5.23 * g / cm3);

  G4Element* elC = nist->FindOrBuildElement("C");
  G4Element* elF = nist->FindOrBuildElement("F");
  auto matPTFE = new G4Material("PTFE", density = 2.2 * g / cm3, nComponents = 2);
  matPTFE->AddElement(elC, nAtoms = 1);
  matPTFE->AddElement(elF, nAtoms = 2);
  auto matFEP = matPTFE;

  G4Material* LN2 = nist->FindOrBuildMaterial("G4_lN2");
  G4Material* LAr = nist->FindOrBuildMaterial("G4_lAr");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //Materials

  G4Material* world_mat = vacuum;
  G4Material* env_mat = LAr;
  G4Material* det_mat;
  if (fMatType == "Enriched") {
	  det_mat = matGeEn;
  }
  else
  {
	  det_mat = matGeNa;
  }

  //     
  // World
  //
  G4double world_size = 20 * cm;
  G4double env_size = 10 * cm;

  G4Box* solidWorld =
	  new G4Box("World",                       //its name
		  0.5 * world_size, 0.5 * world_size, 0.5 * world_size);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_size, 0.5*env_size, 0.5*env_size); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  
  if (fDetGeo == "PPC") {
	  ////////PPC////////
	  //Sphere-like p-point
//     
// Detector
//  
	  double GeChamfer = 2. * mm;
	  double outerGeRadius = 31.35 * mm;
	  double innerGeRadius = 1.50 * mm;
	  double GeHeight1 = 60.80 * mm;
	  double SrcThickness = 0.01 * mm;
	  double lSmallValue = 0.01 * mm;
	  G4double orbradius = 1.89 * mm;

	  // G4Height2 is the depth of the hole for pin contact
	  double GeHeight2 = 4.0 * mm;
	  double GeHeight3 = GeHeight2 + fDetDeadlayer;

	  auto GeP1 = new G4Tubs("GeP1", 0., outerGeRadius, GeHeight1 / 2, 0., twopi);
	  auto GeP2 = new G4Tubs("GeP2", 0., outerGeRadius - GeChamfer, GeChamfer, 0., twopi);
	  auto GeP3 = new G4Torus("GeP3", 0., GeChamfer, outerGeRadius - GeChamfer, 0., twopi);
	  auto GeM = new G4Tubs("solidGeM", 0., innerGeRadius, fDetDeadlayer, 0., twopi);
	  auto GeS = new G4Orb("solidGeS", orbradius);

	  //auto GeCut = new G4Tubs("solidGeCut", 0., 3 * cm, GeHeight2 / 2, 0., twopi);
	  //auto GeM2 = new G4Tubs("solidGeM2", 0., innerGeRadius, GeHeight3 / 2., 0., twopi);

	  G4ThreeVector zTransGe0(0., 0., -GeHeight1 / 2);
	  G4ThreeVector zTransGe1(0., 0., -GeHeight1 / 2 + fDetDeadlayer);
	  G4ThreeVector zTransGe2(0., 0., GeHeight1 / 2 + lSmallValue);

	  // total germanium crystal
	  auto GeTemp1 = new G4UnionSolid("GeTemp1", GeP2, GeP3);
	  auto GeTemp2 = new G4UnionSolid("GeTemp2", GeP1, GeTemp1, 0, zTransGe2);
	  //?
	  auto solidtempTotalCrystal = new G4SubtractionSolid("totaltempCrystal", GeTemp2, GeM, 0, zTransGe0);
	  auto solidTotalCrystal = new G4SubtractionSolid("totalCrystal", solidtempTotalCrystal, GeS, 0, zTransGe1);

	  // bulk
	  auto GeP1In = new G4Tubs("GeP1In", 0., outerGeRadius - fDetDeadlayer, (GeHeight1 - fDetDeadlayer * 2) / 2, 0., twopi);
	  auto GeP2In = new G4Tubs("GeP2In", 0., outerGeRadius - fDetDeadlayer - GeChamfer, GeChamfer, 0., twopi);
	  auto GeP3In = new G4Torus("GeP3In", 0., GeChamfer, outerGeRadius - fDetDeadlayer - GeChamfer, 0., twopi);

	  auto GeMIn = new G4Tubs("solidGeMIn", 0., innerGeRadius, (GeHeight2 - fDetDeadlayer) / 2., 0., twopi);
	  auto GeSIn = new G4Orb("solidGeS", orbradius + fpLayerThickness);

	  auto GeInTemp1 = new G4UnionSolid("GeInTemp1", GeP2In, GeP3In);
	  auto GeInTemp2 = new G4UnionSolid("GeInTemp2", GeP1In, GeInTemp1, 0, G4ThreeVector(0., 0., (GeHeight1 - fDetDeadlayer * 2) / 2.0 + lSmallValue));
	  G4ThreeVector zbulkTrans(0., 0., -(GeHeight1 - fDetDeadlayer * 2 - (GeHeight2 - fDetDeadlayer)) / 2);
	  G4ThreeVector zbulkTransS(0., 0., -GeHeight1 / 2 + fDetDeadlayer);
	  auto bulkCrystal = new G4SubtractionSolid("Bulk", GeInTemp2, GeSIn, 0, zbulkTransS);
	  auto logicBulk = new G4LogicalVolume(bulkCrystal, det_mat, "Bulk");

	  //deadlayer
	  auto solidOuterDeadlayer = new G4SubtractionSolid("OuterDeadlayer", solidTotalCrystal, GeInTemp2, 0, G4ThreeVector(0., 0., 0.));
	  auto logicOuterDeadlayer = new G4LogicalVolume(solidOuterDeadlayer, det_mat, "OuterDeadlayer");

	  //pLayer
	  auto GepLayer = new G4Sphere("solidpLayer", orbradius, orbradius + fpLayerThickness, 0 * pi, 2 * pi, 0 * pi, 0.5 * pi);
	  auto logicGepLayer = new G4LogicalVolume(GepLayer, det_mat, "pLayer");

	  ///////PPC///////
	  new G4PVPlacement(0, G4ThreeVector(), logicBulk, "Bulk", logicEnv, false, 0, checkOverlaps);
	  new G4PVPlacement(0, G4ThreeVector(), logicOuterDeadlayer, "OuterDeadlayer", logicEnv, false, 0, checkOverlaps);
	  new G4PVPlacement(0, zbulkTransS, logicGepLayer, "pLayer", logicEnv, false, 0, checkOverlaps);

	  fScoringVolume = logicBulk;
	  fBulk = logicBulk;
	  fDeadLayer = logicOuterDeadlayer;
	  fpLayer = logicGepLayer;
  }
  
  if (fDetGeo == "flatBEGe") {
	  ////////flatBEGe////////
	  //     
      // Detector
      //  
	  double GeChamfer = 2. * mm;
	  double outerGeRadius = 31.35 * mm;
	  double innerGeRadius = 1.50 * mm;
	  double GeHeight1 = 60.80 * mm;
	  double SrcThickness = 0.01 * mm;
	  double lSmallValue = 0.01 * mm;
	  double groovedepth = 1.5 * mm;
	  double grooveradius = 9 * mm;
	  double groovethickness = 0.001 * mm;
	  double outerplayerthickness = fpLayerThickness;
	  G4double orbradius = 1.89 * mm;

	  // G4Height2 is the depth of the hole for pin contact
	  double GeHeight2 = 4.0 * mm;
	  double GeHeight3 = GeHeight2 + fDetDeadlayer;

	  auto GeP1 = new G4Tubs("GeP1", 0., outerGeRadius, GeHeight1 / 2, 0., twopi);
	  auto GeP2 = new G4Tubs("GeP2", 0., outerGeRadius - GeChamfer, GeChamfer, 0., twopi);
	  auto GeP3 = new G4Torus("GeP3", 0., GeChamfer, outerGeRadius - GeChamfer, 0., twopi);
	  auto GeGroove = new G4Torus("solidgroove", 0., groovedepth, grooveradius, 0., twopi);

	  //auto GeCut = new G4Tubs("solidGeCut", 0., 3 * cm, GeHeight2 / 2, 0., twopi);
	  //auto GeM2 = new G4Tubs("solidGeM2", 0., innerGeRadius, GeHeight3 / 2., 0., twopi);

	  G4ThreeVector zTransGe0(0., 0., -GeHeight1 / 2);
	  G4ThreeVector zTransGe1(0., 0., -GeHeight1 / 2 + fDetDeadlayer);
	  G4ThreeVector zTransGe2(0., 0., GeHeight1 / 2 + lSmallValue);
	  G4ThreeVector zTransGroove(0., 0., -GeHeight1 / 2);

	  // total germanium crystal
	  auto GeTemp1 = new G4UnionSolid("GeTemp1", GeP2, GeP3);
	  auto GeTemp2 = new G4UnionSolid("GeTemp2", GeP1, GeTemp1, 0, zTransGe2);
	  //?
	  //auto solidtempTotalCrystal = new G4SubtractionSolid("totaltempCrystal", GeTemp2, GeGroove, 0, zTransGroove);
	  auto solidTotalCrystal = new G4SubtractionSolid("totalCrystal", GeTemp2, GeGroove, 0, zTransGroove);

	  // bulk
	  auto GeP1In = new G4Tubs("GeP1In", 0., outerGeRadius - fDetDeadlayer, (GeHeight1 - fDetDeadlayer * 2) / 2, 0., twopi);
	  auto GeP2In = new G4Tubs("GeP2In", 0., outerGeRadius - fDetDeadlayer - GeChamfer, GeChamfer, 0., twopi);
	  auto GeP3In = new G4Torus("GeP3In", 0., GeChamfer, outerGeRadius - fDetDeadlayer - GeChamfer, 0., twopi);
	  auto GeP4In = new G4Tubs("GeP4In", 0., grooveradius, fDetDeadlayer / 2, 0., twopi);

	  auto GeLargerGroove = new G4Torus("solidlargergroove", 0., groovedepth + groovethickness, grooveradius, 0., twopi);
	  auto GeOuterpLayer = new G4Tubs("solidouterplayer", 0., grooveradius - groovedepth, outerplayerthickness / 2, 0., twopi);
	  //auto GeOuterpLayer = new G4Tubs("solidouterplayer", 0., grooveradius - groovedepth - groovethickness, outerplayerthickness, 0., twopi);

	  auto GeInTemp1 = new G4UnionSolid("GeInTemp1", GeP2In, GeP3In);
	  auto GeInTemp2 = new G4UnionSolid("GeInTemp2", GeP1In, GeInTemp1, 0, G4ThreeVector(0., 0., (GeHeight1 - fDetDeadlayer * 2) / 2.0 + lSmallValue));
	  auto GeInTemp3 = new G4UnionSolid("GeInTemp3", GeInTemp2, GeP4In, 0, G4ThreeVector(0., 0., -(GeHeight1 - fDetDeadlayer) / 2));
	  G4ThreeVector zbulkTrans(0., 0., -(GeHeight1 - fDetDeadlayer * 2 - (GeHeight2 - fDetDeadlayer)) / 2);
	  G4ThreeVector zbulkTransOuterp(0., 0., -(GeHeight1 - outerplayerthickness) / 2);

	  auto GeInTemp4 = new G4SubtractionSolid("GeInTemp4", GeInTemp3, GeLargerGroove, 0, zTransGroove);
	  auto bulk = new G4SubtractionSolid("Bulk", GeInTemp4, GeOuterpLayer, 0, zbulkTransOuterp);
	  auto logicBulk = new G4LogicalVolume(bulk, det_mat, "Bulk");
	  //auto logicBulk = new G4LogicalVolume(GeInTemp4, det_mat, "Bulk");

	  //deadlayer
	  auto tempdeadlayer=new G4SubtractionSolid("tempdeadlayer", solidTotalCrystal, GeInTemp3, 0, G4ThreeVector(0., 0., 0.));
	  auto solidOuterDeadlayer = new G4SubtractionSolid("OuterDeadlayer", tempdeadlayer, GeLargerGroove, 0, zTransGroove);
	  auto logicOuterDeadlayer = new G4LogicalVolume(solidOuterDeadlayer, det_mat, "OuterDeadlayer");

	  //groove Layer
	  auto GrooveTorus = new G4Torus("solidpLayer", groovedepth, groovedepth + groovethickness, grooveradius, 0., twopi);
	  auto GrooveCut1 = new G4Tubs("groovecut1", 0, grooveradius + groovedepth + groovethickness, groovedepth + groovethickness, 0, twopi);
	  auto GrooveCut2 = new G4Tubs("groovecut2", 0, grooveradius - groovedepth, outerplayerthickness, 0, twopi);
	  auto tempGroove1 = new G4SubtractionSolid("tempgroove1", GrooveTorus, GrooveCut1, 0, G4ThreeVector(0., 0., -(groovedepth + groovethickness)));
	  auto GrooveLayer = new G4SubtractionSolid("GrooveLayer", tempGroove1, GrooveCut2, 0, G4ThreeVector(0., 0., 0.));
	  auto logicGrooveLayer = new G4LogicalVolume(GrooveLayer, det_mat, "GrooveLayer");

	  //outer pLayer
	  auto OuterpLayer = new G4Tubs("solidOuterpLayer", 0, grooveradius - groovedepth, outerplayerthickness / 2, 0, twopi);
	  auto logicOuterpLayer = new G4LogicalVolume(OuterpLayer, det_mat, "OuterpLayer");

	  ///////PPC///////
	  new G4PVPlacement(0, G4ThreeVector(), logicBulk, "Bulk", logicEnv, false, 0, checkOverlaps);
	  new G4PVPlacement(0, G4ThreeVector(), logicOuterDeadlayer, "OuterDeadlayer", logicEnv, false, 0, checkOverlaps);
	  new G4PVPlacement(0, zbulkTransOuterp, logicOuterpLayer, "OuterpLayer", logicEnv, false, 0, checkOverlaps);
	  new G4PVPlacement(0, zTransGroove, logicGrooveLayer, "GrooveLayer", logicEnv, false, 0, checkOverlaps);

	  fScoringVolume = logicBulk;
	  fBulk = logicBulk;
	  fDeadLayer = logicOuterDeadlayer;
	  fpLayer = logicOuterpLayer;
	  
  }



  fEnv = logicEnv;


  ///////Src///////
  //auto GeSrcP1 = new G4Tubs("GeSrcP1", 0., outerGeRadius + SrcThickness, GeHeight1 / 2 + SrcThickness, 0., twopi);
  //auto GeSrcP2 = new G4Tubs("GeSrcP2", 0., outerGeRadius - GeChamfer - SrcThickness, GeChamfer, 0., twopi);
  //auto GeSrcP3 = new G4Torus("GeSrcP3", 0., GeChamfer, outerGeRadius - GeChamfer - SrcThickness, 0., twopi);

  //auto GeSrcTemp1 = new G4UnionSolid("GeSrcTemp1", GeSrcP2, GeSrcP3);
  //G4ThreeVector zTransGeSrc(0., 0., GeHeight1 / 2 + SrcThickness);
  //auto GeSrcTemp2 = new G4UnionSolid("GeSrcTemp2", GeSrcP1, GeSrcTemp1, 0, zTransGeSrc);
  //auto solidSrc = new G4SubtractionSolid("Src", GeSrcTemp2, GeTemp2, 0, G4ThreeVector(0., 0., 0.));
  //auto logicSrc = new G4LogicalVolume(solidSrc, vacuum, "Src");
  //new G4PVPlacement(0, G4ThreeVector(), logicSrc, "Src", logicEnv, false, 0, checkOverlaps);
  








  //Ge p+ surface
  /*
  G4double GeHeight4 = 2. * mm;
  G4ThreeVector zTransGe3(0., 0., -(GeHeight1 - fDetDeadlayer - GeHeight4 - 0.25 * mm) / 2);
  G4double thick_psurface = 1 * um;
  auto SolidGepsurface = new G4Tubs("SolidGepsurface", innerGeRadius - thick_psurface, innerGeRadius, (GeHeight4 - fDetDeadlayer + 0.5 * mm) / 2, 0., twopi);
  auto logicGepsurface = new G4LogicalVolume(SolidGepsurface, LN2, "logicGepsurface");
  new G4PVPlacement(0, zTransGe3, logicGepsurface, "physGepsurface", logicLN2, false, 0, checkOverlaps);
  */
 

  //fSrc = logicSrc;
  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
