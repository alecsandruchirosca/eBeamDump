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
// Authors: Susanna Guatelli, susanna@uow.edu.au,
// Authors: Jeremy Davis, jad028@uowmail.edu.au
//

#include "DetectorConstruction.hh"
#include "globals.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "SensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4UserLimits.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4SystemOfUnits.hh"
#include "Geant4/G4Cons.hh"
#include "Geant4/G4Orb.hh"
#include "Geant4/G4Sphere.hh"
#include "Geant4/G4Trd.hh"
#include "Geant4/G4NistManager.hh"
#include "Geant4/G4Polycone.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSFlatSurfaceFlux.hh"
#include "G4SDParticleFilter.hh"
#include "G4SDChargedFilter.hh"
#include "G4SubtractionSolid.hh"

#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4MagneticField.hh"
#include "G4RotationMatrix.hh"

DetectorConstruction::DetectorConstruction(AnalysisManager* analysis_manager)
{
analysis = analysis_manager;
}

DetectorConstruction::~DetectorConstruction(){

}

G4VPhysicalVolume* DetectorConstruction::Construct2() 
{
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* w_mat = nist->FindOrBuildMaterial("G4_W");
//  G4Material* Ag_mat = nist->FindOrBuildMaterial("G4_Ag");
//  G4Material* target_mat = new G4Material("HDAg", 5*g/cm3, Ag_mat);
  G4Material* C_mat = nist->FindOrBuildMaterial("G4_C");
  G4Material* H_mat = nist->FindOrBuildMaterial("G4_H");
  G4Material* O_mat = nist->FindOrBuildMaterial("G4_O");
  G4Material* target_mat = new G4Material("target", 2.7*g/cm3,4);
  target_mat->AddMaterial(C_mat, 1*perCent);
  target_mat->AddMaterial(H_mat, 4*perCent);
  target_mat->AddMaterial(O_mat, 2*perCent);
  target_mat->AddMaterial(w_mat, 93*perCent);

  G4Material* water = nist->FindOrBuildMaterial("G4_WATER");

  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* vacuum = new G4Material("CETAL", 3.E-9*pascal, air);
  G4Material* glass = nist->FindOrBuildMaterial("G4_GLASS_PLATE");

  G4double density;
  std::vector<G4int> natoms;
  std::vector<G4String> elements;
  

  elements.push_back("C");
  natoms.push_back(5);
  elements.push_back("H");
  natoms.push_back(8);
  elements.push_back("O");
  natoms.push_back(2);
  
  density = 1.190*g/cm3;
   
  G4Material* PMMA = nist->ConstructNewMaterial("PMMA", elements, natoms, density);
  
  elements.clear();
  natoms.clear();
 

  // Construim mediul de propagare
  G4VSolid* worldSolid = new G4Box("world", 300*mm, 300*mm, 2*m);
  // worldLV->SetVisAttributes(G4VisAttributes::Invisible);
  G4VSolid* window = new G4Box("window", 300*mm, 300*mm, 10*mm);
  G4VSolid* airSolid = new G4Box("Air", 300*mm, 300*mm, 0.5*m);
  G4VSolid* layer1 = new G4Box("L1", 300*mm, 300*mm, 5*cm);
  G4VSolid* layer2 = new G4Box("L2", 300*mm, 300*mm, 2.5*cm);
  G4VSolid* layer3 = new G4Box("L3", 300*mm, 300*mm, 5*cm);
  G4VSolid* layer4 = new G4Box("L4", 300*mm, 300*mm, 2.5*cm);
// Logical Volumes
  G4LogicalVolume* worldLogical = new G4LogicalVolume(worldSolid, vacuum, "WorldLV");
  G4LogicalVolume* windowLV = new G4LogicalVolume(window, glass, "WindowLV");
  G4LogicalVolume* airLV = new G4LogicalVolume(airSolid, air, "AirSolid");
  G4LogicalVolume* layer1LV = new G4LogicalVolume(layer1, PMMA, "L1PV");
  G4LogicalVolume* layer2LV = new G4LogicalVolume(layer2, PMMA, "L2PV");
  G4LogicalVolume* layer3LV = new G4LogicalVolume(layer3, PMMA, "L3PV");
  G4LogicalVolume* layer4LV = new G4LogicalVolume(layer4, PMMA, "L41PV");

  G4VisAttributes* preAttr = new G4VisAttributes(G4Color(0.,0.2,0.));
  preAttr->SetForceWireframe(true);
  preAttr->SetForceAuxEdgeVisible(true);
  // layer1LV->SetVisAttributes(preAttr);

  G4VisAttributes* cilAttr = new G4VisAttributes(G4Color(1.,0.2,0.));
  cilAttr->SetForceWireframe(true);
  cilAttr->SetForceAuxEdgeVisible(true);


  // detLV->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* detAttr = new G4VisAttributes(G4Color(1.,0.2,0.));
  detAttr->SetForceWireframe(true);
  detAttr->SetForceAuxEdgeVisible(true);
  // airLV->SetVisAttributes(preAttr);


  G4VisAttributes* trgAttr = new G4VisAttributes(G4Color(0., 1., 0.));
  trgAttr->SetForceAuxEdgeVisible(true);
  trgAttr->SetForceWireframe(true);
  // targetLV->SetVisAttributes(trgAttr);

  G4VisAttributes* foilAttr = new G4VisAttributes(G4Color(1., 0., 0.));
  foilAttr->SetForceAuxEdgeVisible(true);
  foilAttr->SetForceWireframe(true);
  layer4LV->SetVisAttributes(foilAttr);

  // foilLV->SetVisAttributes(foilAttr);

  G4bool debugGeom = true;

// Placements
  G4VPhysicalVolume* worldPhysical = new G4PVPlacement(0, G4ThreeVector(), worldLogical, "WorldPhysical", 0, false, 0, false);
  new G4PVPlacement(0, G4ThreeVector(0, 0*mm, 0.4*m), windowLV, "WindowPhys", worldLogical, false, 0, debugGeom);
  new G4PVPlacement(0, G4ThreeVector(0*mm, 0*mm, 0.9*m), airLV, "AirPhys", worldLogical, false, 0, debugGeom);
  new G4PVPlacement(0, G4ThreeVector(0*mm, 0*mm, 1.45*m), layer1LV, "Layer1Phys", worldLogical, false, 0, debugGeom);
  new G4PVPlacement(0, G4ThreeVector(0*mm, 0*mm, 1.525*m), layer2LV, "Layer2Phys", worldLogical, false, 0, debugGeom);
  new G4PVPlacement(0, G4ThreeVector(0*mm, 0*mm, 1.6*m), layer3LV, "Layer3Phys", worldLogical, false, 0, debugGeom);
  new G4PVPlacement(0, G4ThreeVector(0*mm, 0*mm, 1.675*m), layer4LV, "Layer4Phys", worldLogical, false, 0, debugGeom);

//  worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
 
  return worldPhysical;
}



G4VPhysicalVolume* DetectorConstruction::Construct() 
{
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* w_mat = nist->FindOrBuildMaterial("G4_W");
//  G4Material* Ag_mat = nist->FindOrBuildMaterial("G4_Ag");
//  G4Material* target_mat = new G4Material("HDAg", 5*g/cm3, Ag_mat);
  G4Material* C_mat = nist->FindOrBuildMaterial("G4_C");
  G4Material* H_mat = nist->FindOrBuildMaterial("G4_H");
  G4Material* O_mat = nist->FindOrBuildMaterial("G4_O");
  G4Material* target_mat = new G4Material("target", 2.7*g/cm3,4);
  target_mat->AddMaterial(C_mat, 1*perCent);
  target_mat->AddMaterial(H_mat, 4*perCent);
  target_mat->AddMaterial(O_mat, 2*perCent);
  target_mat->AddMaterial(w_mat, 93*perCent);

  G4Material* water = nist->FindOrBuildMaterial("G4_WATER");

  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* vacuum = new G4Material("CETAL", 3.E-9*pascal, air);
  G4Material* glass = nist->FindOrBuildMaterial("G4_GLASS_PLATE");

  G4double density;
  std::vector<G4int> natoms;
  std::vector<G4String> elements;
  

  elements.push_back("C");
  natoms.push_back(5);
  elements.push_back("H");
  natoms.push_back(8);
  elements.push_back("O");
  natoms.push_back(2);
  
  density = 1.190*g/cm3;
   
  G4Material* PMMA = nist->ConstructNewMaterial("PMMA", elements, natoms, density);
  
  elements.clear();
  natoms.clear();
 

  // Construim mediul de propagare
  G4VSolid* worldSolid = new G4Box("world", 300*mm, 300*mm, 2.5*m);
  // worldLV->SetVisAttributes(G4VisAttributes::Invisible);
  G4VSolid* window = new G4Box("window", 300*mm, 300*mm, 10*mm);
  G4VSolid* airSolid = new G4Box("Air", 300*mm, 300*mm, 0.5*m);
  G4VSolid* layer1 = new G4Box("L1", 300*mm, 300*mm, 15*cm);
  G4VSolid* layer2 = new G4Box("L2", 300*mm, 300*mm, 2.5*cm);
  G4VSolid* layer3 = new G4Box("L3", 300*mm, 300*mm, 15*cm);
  G4VSolid* layer4 = new G4Box("L4", 300*mm, 300*mm, 2.5*cm);
// Logical Volumes
  G4LogicalVolume* worldLogical = new G4LogicalVolume(worldSolid, vacuum, "WorldLV");
  G4LogicalVolume* windowLV = new G4LogicalVolume(window, glass, "WindowLV");
  G4LogicalVolume* airLV = new G4LogicalVolume(airSolid, air, "AirSolid");
  G4LogicalVolume* layer1LV = new G4LogicalVolume(layer1, water, "L1PV");
  G4LogicalVolume* layer2LV = new G4LogicalVolume(layer2, water, "L2PV");
  G4LogicalVolume* layer3LV = new G4LogicalVolume(layer3, water, "L3PV");
  G4LogicalVolume* layer4LV = new G4LogicalVolume(layer4, water, "L41PV");
  worldLogical->SetVisAttributes(G4VisAttributes::Invisible);

  G4bool debugGeom = true;
// Placements
  G4VPhysicalVolume* worldPhysical = new G4PVPlacement(0, G4ThreeVector(), worldLogical, "WorldPhysical", 0, false, 0, false);
  new G4PVPlacement(0, G4ThreeVector(0, 0*mm, 0.4*m), windowLV, "WindowPhys", worldLogical, false, 0, debugGeom);
  new G4PVPlacement(0, G4ThreeVector(0*mm, 0*mm, 0.9*m), airLV, "AirPhys", worldLogical, false, 0, debugGeom);
  new G4PVPlacement(0, G4ThreeVector(0*mm, 0*mm, 1.4*m+15*cm), layer1LV, "Layer1Phys", worldLogical, false, 0, debugGeom);
  new G4PVPlacement(0, G4ThreeVector(0*mm, 0*mm, 1.7*m+2.5*cm), layer2LV, "Layer2Phys", worldLogical, false, 0, debugGeom);
  new G4PVPlacement(0, G4ThreeVector(0*mm, 0*mm, 1.75*m+15*cm), layer3LV, "Layer3Phys", worldLogical, false, 0, debugGeom);
  new G4PVPlacement(0, G4ThreeVector(0*mm, 0*mm, 2.05*m+2.5*cm), layer4LV, "Layer4Phys", worldLogical, false, 0, debugGeom);

//  worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
 
  return worldPhysical;
}


void DetectorConstruction::ConstructSDandField()
{
  // G4FieldManager* fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  // G4MagneticField* magField= new G4UniformMagField(G4ThreeVector(0.6*tesla,0.,0.));
  // fieldManager->SetDetectorField(magField);
  // fieldManager->CreateChordFinder(magField);
  
  // G4FieldManager* globalFieldMgr= G4TransportationManager::GetTransportationManager()->GetFieldManager();
  // globalFieldMgr->SetDetectorField(magField);

  /*
   G4MultiFunctionalDetector* myScorer = new G4MultiFunctionalDetector("ProtonScorer");
   G4VPrimitiveScorer* totalFlux = new G4PSFlatSurfaceFlux("TotalSurfaceFlux",0);
   myScorer->RegisterPrimitive(totalFlux);
   G4VPrimitiveScorer* protonSurfFlux = new G4PSFlatSurfaceFlux("ProtonSurfaceFlux", 0);
   G4SDParticleFilter* protonFilter = new G4SDParticleFilter("protonFilter");
   protonFilter->add("proton");
   protonSurfFlux->SetFilter(protonFilter);
   myScorer->RegisterPrimitive(protonSurfFlux);
*/
}
