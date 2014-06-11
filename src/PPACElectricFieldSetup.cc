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
// $Id: PPACElectricFieldSetup.cc,v 1.4 2008-05-14 15:27:29 tnikitin Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//   User Field class implementation.
//

#include "PPACElectricFieldSetup.hh"
#include "PPACFieldMessenger.hh"

#include "G4UniformElectricField.hh"
#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"

#include "G4RKG3_Stepper.hh"

//////////////////////////////////////////////////////////////////////////
//
//  Constructors:

PPACElectricFieldSetup::PPACElectricFieldSetup()
  : fChordFinder(0), fLocalChordFinder(0), fStepper(0), fIntgrDriver(0)
{
  fEMfield = new G4UniformElectricField(
                 G4ThreeVector(0.0*kilovolt/mm,0.0*kilovolt/mm,0.0*volt/mm));

  fLocalEMfield = new G4UniformElectricField(
                 G4ThreeVector(0.0*kilovolt/mm,0.0*kilovolt/mm,0.0*volt/mm));

  fFieldMessenger = new PPACFieldMessenger(this) ;  

  fEquation = new G4EqMagElectricField(fEMfield); 
  fLocalEquation = new G4EqMagElectricField(fLocalEMfield);

  fMinStep     = 0.010*mm ; // minimal step of 10 microns
  fLocalMinStep     = 0.010*mm ; // minimal step of 10 microns
  fStepperType = 4 ;        // ClassicalRK4 -- the default stepper

  fFieldManager = GetGlobalFieldManager();
  fLocalFieldManager = new G4FieldManager();

  UpdateField();

}

/////////////////////////////////////////////////////////////////////////////////

PPACElectricFieldSetup::PPACElectricFieldSetup(G4ThreeVector pFV)
  : fChordFinder(0), 
    fLocalChordFinder(0),
    fStepper(0),
    fIntgrDriver(0)
{    
  fEMfield = new G4UniformElectricField(pFV);
//  GetGlobalFieldManager()->CreateChordFinder(fEMfield);
  fLocalEMfield = new G4UniformElectricField(pFV);

  fFieldMessenger = new PPACFieldMessenger(this) ; 
 
  fEquation = new G4EqMagElectricField(fEMfield); 
  fLocalEquation = new G4EqMagElectricField(fLocalEMfield);
  fMinStep     = 0.010*mm ; // minimal step of 10 microns
  fLocalMinStep     = 0.010*mm ; // minimal step of 10 microns
  fStepperType = 4 ;        // ClassicalRK4 -- the default stepper

  fFieldManager = GetGlobalFieldManager();
  fLocalFieldManager = new G4FieldManager();

  UpdateField();
}


////////////////////////////////////////////////////////////////////////////////

PPACElectricFieldSetup::~PPACElectricFieldSetup()
{
  if(fChordFinder) delete fChordFinder;
  if(fStepper)     delete fStepper;
  if(fEquation)    delete fEquation;   
  if(fEMfield)     delete fEMfield;
}

/////////////////////////////////////////////////////////////////////////////
//
// Register this field to 'global' Field Manager and 
// Create Stepper and Chord Finder with predefined type, minstep (resp.)
//

void PPACElectricFieldSetup::UpdateField()
{
  SetStepper();

  G4cout<<"The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

  fFieldManager->SetDetectorField(fEMfield );
  fLocalFieldManager->SetDetectorField(fLocalEMfield );
//  G4Box* EField_box = new G4Box("EField_box", 100.*mm/2, 100.*mm/2, 4.*mm/2);
//  G4LogicalVolume* EField_log = new G4LogicalVolume(EField_box, Universal_Vacuum, "EField_log", 0, 0, 0);

  if(fChordFinder) delete fChordFinder;
  if(fLocalChordFinder) delete fLocalChordFinder;

  // fChordFinder = new G4ChordFinder( fEMfield, fMinStep, fStepper);
  fIntgrDriver = new G4MagInt_Driver(fMinStep, 
				     fStepper, 
                                     fStepper->GetNumberOfVariables() );
  G4cout<<"test2 "<<G4endl ;

  fLocalIntgrDriver = new G4MagInt_Driver(fMinStep,
                                     fStepper,
                                     fStepper->GetNumberOfVariables() );
  G4cout<<"test3 "<<G4endl ;
  fChordFinder = new G4ChordFinder(fIntgrDriver);
  fLocalChordFinder = new G4ChordFinder(fLocalIntgrDriver);

  fFieldManager->SetChordFinder( fChordFinder );
  fLocalFieldManager->SetChordFinder( fLocalChordFinder );
}

/////////////////////////////////////////////////////////////////////////////
//
// Set stepper according to the stepper type
//

void PPACElectricFieldSetup::SetStepper()
{
  G4int nvar = 8;

  if(fStepper) delete fStepper;
  
  switch ( fStepperType ) 
  {
    case 0:  
      fStepper = new G4ExplicitEuler( fEquation, nvar );
      fLocalStepper = new G4ExplicitEuler( fLocalEquation, nvar ); 
      G4cout<<"G4ExplicitEuler is calledS"<<G4endl;     
      break;
    case 1:  
      fStepper = new G4ImplicitEuler( fEquation, nvar );
      fLocalStepper = new G4ImplicitEuler( fLocalEquation, nvar );      
      G4cout<<"G4ImplicitEuler is called"<<G4endl;     
      break;
    case 2:  
      fStepper = new G4SimpleRunge( fEquation, nvar );        
      fLocalStepper = new G4SimpleRunge( fLocalEquation, nvar );
      G4cout<<"G4SimpleRunge is called"<<G4endl;     
      break;
    case 3:  
      fStepper = new G4SimpleHeum( fEquation, nvar );
      fLocalStepper = new G4SimpleHeum( fLocalEquation, nvar );         
      G4cout<<"G4SimpleHeum is called"<<G4endl;     
      break;
    case 4:
      fStepper = new G4ClassicalRK4( fEquation, nvar );
      G4cout<<"G4ClassicalRK4 (default) is called"<<G4endl;
      break;
    case 5:
//      fStepper = new G4CashKarpRKF45( fEquation, nvar );
//      G4cout<<"G4CashKarpRKF45 is called"<<G4endl;
      fStepper = 0;
      break;
    case 6:  
      fStepper = 0; // new G4RKG3_Stepper( fEquation, nvar );       
      G4cout<<"G4RKG3_Stepper is not currently working for Electric Field"<<G4endl;     
      break;
    case 7:  
      fStepper = 0; // new G4HelixExplicitEuler( fEquation ); 
      G4cout<<"G4HelixExplicitEuler is not valid for Electric Field"<<G4endl;     
      break;
    case 8:  
      fStepper = 0; // new G4HelixImplicitEuler( fEquation ); 
      G4cout<<"G4HelixImplicitEuler is not valid for Electric Field"<<G4endl;     
      break;
    case 9:  
      fStepper = 0; // new G4HelixSimpleRunge( fEquation );   
      G4cout<<"G4HelixSimpleRunge is not valid for Electric Field"<<G4endl;     
      break;
    default: fStepper = 0;
  }
}

/////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field to fieldValue along Z
//

void PPACElectricFieldSetup::SetFieldValue(G4double fieldValue)
{
  G4ThreeVector fieldVector( 0.0, 0.0, fieldValue );  

  SetFieldValue( fieldVector );
}

///////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field value to fieldVector
//

void PPACElectricFieldSetup::SetFieldValue(G4ThreeVector fieldVector)
{
  // Find the Field Manager for the global field
  G4FieldManager* fieldMgr= GetGlobalFieldManager();
    
  if(fieldVector != G4ThreeVector(0.,0.,0.))
  { 
    if(fEMfield) delete fEMfield;
    fEMfield = new  G4UniformElectricField(fieldVector);

    fEquation->SetFieldObj(fEMfield);  // must now point to the new field

    // UpdateField();
   
    fieldMgr->SetDetectorField(fEMfield);
//    gapLV->SetFieldManager(fieldMgr,true);
  }
  else 
  {
    // If the new field's value is Zero, then it is best to
    //  insure that it is not used for propagation.
    if(fEMfield) delete fEMfield;
    fEMfield = 0;
    fEquation->SetFieldObj(fEMfield);   // As a double check ...

    G4MagneticField* fEMfield = 0;
     fieldMgr->SetDetectorField(fEMfield);
  }
}

////////////////////////////////////////////////////////////////////////////////
//
//  Utility method

G4FieldManager*  PPACElectricFieldSetup::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()
	 ->GetFieldManager();
}
