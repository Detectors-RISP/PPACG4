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
// $Id$
// 
/// \file PPACDetectorConstruction.cc
/// \brief Implementation of the PPACDetectorConstruction class

#include "PPACDetectorConstruction.hh"
#include "PPACDetectorMessenger.hh"
#include "PPACCalorimeterSD.hh"
#include "PPACElectricFieldSetup.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformElectricField.hh"
#include "G4UniformMagField.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include <stdio.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PPACDetectorConstruction::PPACDetectorConstruction()
 : G4VUserDetectorConstruction(),
   fEmFieldSetup(0),
   fMagField(0),
   fCheckOverlaps(true)
{
  fMessenger = new PPACDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PPACDetectorConstruction::~PPACDetectorConstruction()
{ 
  delete fMagField;
  if (fEmFieldSetup) delete fEmFieldSetup;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* PPACDetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();

  fEmFieldSetup = new PPACElectricFieldSetup();
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PPACDetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;
  nistManager->FindOrBuildMaterial("G4_Pb", fromIsotopes);
  
  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density, massfraction; 
//  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, STP_Temperature+15.*kelvin, 3.e-18*pascal);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  G4double A, Z;
  G4int nel, natoms;
  G4String name, symbol;
//  const double inch = 2.54*cm;
//  char name1[30],name2[30];

  G4Colour aqua(0.247, 0.8, 1.0);
  G4Colour magenda(1.0, 0.0, 1.0);
  G4Colour yellow(1.0, 1.0, 0.0);
  G4Colour red(1.0, 0.0, 0.0);

  const G4double expTemp= STP_Temperature+15.*kelvin; // temperature of experimetal hall is 15C
  const G4double isobutanedatatmp= STP_Temperature+15.*kelvin; //2.51 mg mL-` , at 15C, 100kPa
  const G4double Perfluoropropanedatatmp = STP_Temperature;

  const G4double pressure = 1.3158e-2*atmosphere; //6 torr = 0.007895 atm, 10 torr = 0.013158, 50 = 0.065789
  const G4double normalatm = 1*atmosphere;
  A= 1.00794 *g/mole;
  G4Element* elH= new G4Element(name="Hydrogen", symbol="H", Z=1., A);

  A= 18.9984032 *g/mole;
  G4Element* elF= new G4Element(name="Fluorine", symbol="F", Z=9., A);

  A= 12.011 *g/mole;
  G4Element* elC= new G4Element(name="Carbon", symbol="C", Z=6., A);

  A= 15.9994 *g/mole;
  G4Element* elO= new G4Element(name="Oxygen", symbol="O", Z=8., A); 

  A= 107.8682 *g/mole; // Genie
  G4Element* elAg= new G4Element(name="Agsilver", symbol="O", Z=47., A);

  A= 26.9815386 *g/mole;
  G4Element* elAl = new G4Element(name="Aluminum", symbol="Al", Z=13., A);

//  const G4double denIsobutene = 1.356e-3*g/cm3 * STP_Temperature/expTemp;
//  const G4double denIsobutene = 2.51e-3*g/cm3 * STP_Temperature/expTemp; //2.51 mg mL-` , at 15C, 100kPa
  const G4double denIsobutene = 2.51e-3*g/cm3 * isobutanedatatmp/expTemp*pressure/normalatm; 

  G4Material* Isobutene = new G4Material(name="Isobutene", denIsobutene, nel = 2, kStateGas, expTemp, pressure);
  Isobutene->AddElement(elC, natoms=4);
  Isobutene->AddElement(elH, natoms=10);

  const G4double denPerfluoropropane = 8.17e-3*g/cm3 * Perfluoropropanedatatmp/expTemp*pressure/normalatm;
  G4Material* Perfluoropropane = new G4Material(name="Perfluoropropane", denPerfluoropropane, nel = 2, kStateGas, expTemp, pressure);
  Perfluoropropane->AddElement(elC, natoms=3);
  Perfluoropropane->AddElement(elF, natoms=8);

  G4cout << denIsobutene;
//  G4cout << *(G4material::GetMaterialTable());

   // Mylar as Polyethylene terephthalate
  density = 1.38 *g/cm3;
  G4Material* Mylar = new G4Material("Mylar", density, nel = 3);
  Mylar -> AddElement(elC, 10);
  Mylar -> AddElement(elH, 8);
  Mylar -> AddElement(elO, 4);

  density = 10.49 *g/cm3;
  G4Material* Silver = new G4Material("Silver", density, nel = 1);
  Silver -> AddElement(elAg, massfraction=1.0);

  density = 2.70 *g/cm3;
  G4Material* Aluminium = new G4Material("Aluminium", density, nel = 1);
  Aluminium -> AddElement(elAl, massfraction=1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* PPACDetectorConstruction::DefineVolumes()
{
  // Geometry parameters
//  G4int nofLayers = 1;
  G4double mylarThickness = 2.0e-3*mm; //mylar
  G4double gapThickness  = 4.*mm;
  G4double gap2Thickness  = 13.*mm;
  G4double detectorSizeXY   = 10.*cm;
  G4double electrodeThickness   = 3.0e-5*mm; // Ag Silver
  G4double aluminiumThickness   = 2.0e-4*mm; //Aluminium 200nm

//  G4double layerThickness = 3.0 * mylarThickness + 2.0 * gapThickness + 3.0 * electrodeThickness;
  G4double detectorThickness = 3.0 * mylarThickness + 2.0 * gapThickness + 3.0 * electrodeThickness;
  G4double detectorThickness2 = 5.0 * mylarThickness + 2.0 * (gapThickness+gap2Thickness+aluminiumThickness) + 3.0 * electrodeThickness;
  G4double worldSizeXY = 4 * detectorSizeXY;
  G4double worldSizeZ  = 4 * detectorThickness2; 
  
  // Get materials
  G4Material* defaultMaterial = G4Material::GetMaterial("Galactic");
  G4Material* mylarMaterial = G4Material::GetMaterial("Mylar");
  G4Material* electrodeMaterial = G4Material::GetMaterial("Silver");
  G4Material* gapMaterial = G4Material::GetMaterial("Isobutene");
  G4Material* gap2Material = G4Material::GetMaterial("Isobutene");
  G4Material* surfaceMaterial = G4Material::GetMaterial("Aluminium");
  
  if ( ! defaultMaterial || ! mylarMaterial || ! gapMaterial || ! electrodeMaterial || ! gap2Material || ! surfaceMaterial) {
    G4cerr << "Cannot retrieve materials already defined. " << G4endl;
    G4cerr << "Exiting application " << G4endl;
    exit(1);
  }  
  
// Print materials
   G4cout << *(G4Material::GetMaterialTable()) << G4endl;
 
  //     
  // World
  //
  G4VSolid* worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                         
  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                               
  // Detector
  //  
  G4VSolid* detectorS
    = new G4Box("detector",     // its name
                 detectorSizeXY/2, detectorSizeXY/2, detectorThickness2/2); // its size
                         
  G4LogicalVolume* detectorLV
    = new G4LogicalVolume(
                 detectorS,     // its solid
                 defaultMaterial,  // its material
                 "detector");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 detectorLV,          // its logical volume                         
                 "detector",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
   
  //                                 
  // Layer
  //
/*  G4VSolid* layerS 
    = new G4Box("Layer",           // its name
                 detectorSizeXY/2, detectorSizeXY/2, layerThickness/2); //its size
                         
  G4LogicalVolume* layerLV
    = new G4LogicalVolume(
                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Layer");         // its name

  new G4PVReplica(
                 "Layer",          // its name
                 layerLV,          // its logical volume
                 detectorLV,          // its mother
                 kZAxis,           // axis of replication
                 nofLayers,        // number of replica
                 layerThickness);  // witdth of replica
  
*/
  //                               
  // Mylar
  //
  G4VSolid* surfaceS
    = new G4Box("surface",            // its name
                 detectorSizeXY/2, detectorSizeXY/2, aluminiumThickness/2); // its size

  G4LogicalVolume* surfaceLV
    = new G4LogicalVolume(
                 surfaceS,        // its solid
                 surfaceMaterial, // its material
                 "surface");          // its name

   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., aluminiumThickness/2 - detectorThickness2/2), // its position
                 surfaceLV,       // its logical volume
                 "surface",           // its name
                 detectorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., detectorThickness2/2 - aluminiumThickness/2), // its position
                 surfaceLV,       // its logical volume
                 "surface",           // its name
                 detectorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlap


  G4VSolid* mylarS 
    = new G4Box("myl",            // its name
                 detectorSizeXY/2, detectorSizeXY/2, mylarThickness/2); // its size
                         
  G4LogicalVolume* mylarLV
    = new G4LogicalVolume(
                 mylarS,        // its solid
                 mylarMaterial, // its material
                 "myl");          // its name

  G4int nofgap = 2;

   for(G4int i=0; i<nofgap+1; i++)
   {                                   
   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., mylarThickness/2 + electrodeThickness + i*(electrodeThickness+gapThickness+mylarThickness) - detectorThickness/2), // its position
                 mylarLV,       // its logical volume                         
                 "myl",           // its name
                 detectorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  }

  //for surface mylar
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., aluminiumThickness + mylarThickness/2 - detectorThickness2/2), // its position
                 mylarLV,       // its logical volume
                 "myl",           // its name
                 detectorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., detectorThickness2/2 - aluminiumThickness -  mylarThickness/2), // its position
                 mylarLV,       // its logical volume
                 "myl",           // its name
                 detectorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlap

  //
  // electrode
  //
  G4VSolid* electroS
    = new G4Box("electro",            // its name
                 detectorSizeXY/2, detectorSizeXY/2, electrodeThickness/2); // its size

  G4LogicalVolume* electroLV
    = new G4LogicalVolume(
                 electroS,        // its solid
                 electrodeMaterial, // its material
                 "electro");          // its name

   for(G4int i=0; i<nofgap+1; i++)
   {
   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., electrodeThickness/2 + i*(electrodeThickness+gapThickness+mylarThickness) - detectorThickness/2), // its position
                 electroLV,       // its logical volume
                 "electro",           // its name
                 detectorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  }


  //                               
  // Gap
  //
  G4VSolid* gapS 
    = new G4Box("Gap",             // its name
                 detectorSizeXY/2, detectorSizeXY/2, gapThickness/2); // its size
                         
  G4LogicalVolume* gapLV
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "Gap");           // its name

  // Set local field manager and local field in radiator and its daughters:
  //
  G4bool allLocal = true ;
  gapLV->SetFieldManager( fEmFieldSetup->GetLocalFieldManager(),
                                  allLocal ) ;

   for(G4int i=0; i<nofgap; i++)
   {                                   
   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., gapThickness/2 + mylarThickness +electrodeThickness + i*(electrodeThickness+gapThickness+mylarThickness) - detectorThickness/2), // its position
                 gapLV,            // its logical volume                         
                 "Gap",            // its name
                 detectorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
 
  }


   G4VSolid* gap2S
    = new G4Box("Gap2",             // its name
                 detectorSizeXY/2, detectorSizeXY/2, gap2Thickness/2); // its size

  G4LogicalVolume* gap2LV
    = new G4LogicalVolume(
                 gap2S,             // its solid
                 gap2Material,      // its material
                 "Gap2");           // its name

//  G4bool allLocal = true ;
//  gap2LV->SetFieldManager( fEmFieldSetup->GetLocalFieldManager(),
//                                  allLocal ) ;

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., aluminiumThickness + mylarThickness + gap2Thickness/2 - detectorThickness2/2), // its position
                 gap2LV,       // its logical volume
                 "Gap2",           // its name
                 detectorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., detectorThickness2/2 - aluminiumThickness -  mylarThickness - gap2Thickness/2), // its position
                 gap2LV,       // its logical volume
                 "Gap2",           // its name
                 detectorLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlap

  //
  // print parameters
  //
  G4cout << "\n------------------------------------------------------------"
         << "\n---> The calorimeter is " << nofgap << " gaps of: [ "
         << mylarThickness/mm << "mm of " << mylarMaterial->GetName() 
         << " + "
         << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " 
         << "\n------------------------------------------------------------\n";
  
  
  // 
  // Sensitive detectors
  //

  PPACCalorimeterSD* mylarSD 
    = new PPACCalorimeterSD("MylarSD", "MylarHitsCollection", nofgap+1+2);
  G4SDManager::GetSDMpointer()->AddNewDetector(mylarSD );
  mylarLV->SetSensitiveDetector(mylarSD);

  PPACCalorimeterSD* electroSD
    = new PPACCalorimeterSD("ElectroSD", "ElectroHitsCollection", nofgap+1);
  G4SDManager::GetSDMpointer()->AddNewDetector(electroSD );
  electroLV->SetSensitiveDetector(electroSD);

  PPACCalorimeterSD* gapSD 
    = new PPACCalorimeterSD("GapSD", "GapHitsCollection", nofgap);
  G4SDManager::GetSDMpointer()->AddNewDetector(gapSD );
  gapLV->SetSensitiveDetector(gapSD);
  
  PPACCalorimeterSD* gap2SD
    = new PPACCalorimeterSD("Gap2SD", "Gap2HitsCollection", nofgap);
  G4SDManager::GetSDMpointer()->AddNewDetector(gap2SD );
  gap2LV->SetSensitiveDetector(gap2SD);

  PPACCalorimeterSD* surfaceSD
    = new PPACCalorimeterSD("SurfaceSD", "SurfaceHitsCollection", nofgap);
  G4SDManager::GetSDMpointer()->AddNewDetector(surfaceSD );
  surfaceLV->SetSensitiveDetector(surfaceSD); 

  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::Invisible);

  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  detectorLV->SetVisAttributes(simpleBoxVisAtt);

//  EField();

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PPACDetectorConstruction::SetMagField(G4double fieldValue)
{
  // Apply a global uniform magnetic field along X axis
  G4FieldManager* fieldManager
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  // Delete the existing magnetic field
  if ( fMagField )  delete fMagField; 

  if ( fieldValue != 0. ) {
    // create a new one if not null
    fMagField 
      = new G4UniformMagField(G4ThreeVector(fieldValue, 0., 0.));
      
    fieldManager->SetDetectorField(fMagField);
    fieldManager->CreateChordFinder(fMagField);
  } 
  else {
    fMagField = 0;
    fieldManager->SetDetectorField(fMagField);
  }
}

/*
void  PPACDetectorConstruction::EField()
{
  // Apply a global uniform magnetic field along X axis
  G4FieldManager* fieldManager
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  // Delete the existing electric field
  if ( fEField )  delete fEField;

    // create a new one if not null
    fEField
      = new G4UniformElectricField(G4ThreeVector(0., 0., 200.0*volt/mm));
    fEquation->SetFieldObj(fEfield);  // must now point to the new field
    
    fieldManager->SetDetectorField(fEField);
//    fieldManager->CreateChordFinder(fEField);
  }
}

void  PPACDetectorConstruction::SetEField(G4double efieldValue)
{
  // Apply a global uniform magnetic field along X axis
  G4FieldManager* fieldManager
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  // Delete the existing electric field
  if ( fEField )  delete fEField;

  if ( efieldValue != 0. ) {
    // create a new one if not null
    fEField
      = new G4UniformElectricField(G4ThreeVector(0., 0., efieldValue));
    fEquation->SetFieldObj(fEfield);  // must now point to the new field

    fieldManager->SetDetectorField(fEField);
//    fieldManager->CreateChordFinder(fEField);
  }
  else {
    if(fEfield) delete fEfield;
    fEField = 0;
    fEquation->SetFieldObj(fEfield);   // As a double check ...

    G4MagneticField* fEfield = 0;
    fieldManager->SetDetectorField(fEField);
  }
}*/
