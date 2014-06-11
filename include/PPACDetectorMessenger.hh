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
/// \file PPACDetectorMessenger.hh
/// \brief Definition of the PPACDetectorMessenger class

#ifndef PPACDetectorMessenger_h
#define PPACDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PPACDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;

/// Messenger class that defines commands for PPACDetectorConstruction.
///
/// It implements commands:
/// - /B4/det/setMagField value unit

class PPACDetectorMessenger: public G4UImessenger
{
  public:
    PPACDetectorMessenger(PPACDetectorConstruction* detectorConstruction);
    virtual ~PPACDetectorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    PPACDetectorConstruction*  fDetectorConstruction;
    
    G4UIdirectory*             fB4Directory;
    G4UIdirectory*             fDetDirectory;
    G4UIcmdWithADoubleAndUnit* fSetMagFieldCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

