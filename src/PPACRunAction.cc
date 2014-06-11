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
/// \file PPACRunAction.cc
/// \brief Implementation of the PPACRunAction class

#include "PPACRunAction.hh"
#include "PPACAnalysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PPACRunAction::PPACRunAction()
 : G4UserRunAction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PPACRunAction::~PPACRunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PPACRunAction::BeginOfRunAction(const G4Run* run)
{ 
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;

  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Book histograms, ntuple
  //
  
  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in PPACAnalysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  // Open an output file
  //
  G4String fileName = "PPAC";
  analysisManager->OpenFile(fileName);
  analysisManager->SetFirstHistoId(1);

  // Creating histograms
  //
  analysisManager->CreateH1("1","Edep in mylar", 100, 0., 800*MeV);
  analysisManager->CreateH1("2","Edep in gap", 100, 0., 100*MeV);
  analysisManager->CreateH1("3","Edep in electrode", 100, 0., 100*MeV);
  analysisManager->CreateH1("4","Edep in gap2", 100, 0., 100*MeV);
  analysisManager->CreateH1("5","Edep in surface", 100, 0., 100*MeV);
  analysisManager->CreateH1("6","trackL in mylar", 100, 0., 50*cm);
  analysisManager->CreateH1("7","trackL in gap", 100, 0., 50*cm);
  analysisManager->CreateH1("8","trackL in electrode", 100, 0., 50*cm);
  analysisManager->CreateH1("9","trackL in gap2", 100, 0., 50*cm);
  analysisManager->CreateH1("10","trackL in surface", 100, 0., 50*cm);
 
  // Creating ntuple
  //
  analysisManager->CreateNtuple("B6", "Edep and TrackL");
  analysisManager->CreateNtupleDColumn("Emylar");
  analysisManager->CreateNtupleDColumn("Egap");
  analysisManager->CreateNtupleDColumn("Eelectrode");
  analysisManager->CreateNtupleDColumn("Egap2");
  analysisManager->CreateNtupleDColumn("Esurface");
  analysisManager->CreateNtupleDColumn("Lmylar");
  analysisManager->CreateNtupleDColumn("Lgap");
  analysisManager->CreateNtupleDColumn("Lelectrode");
  analysisManager->CreateNtupleDColumn("Lgap2");
  analysisManager->CreateNtupleDColumn("Lsurface");
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PPACRunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int nofEvents = aRun->GetNumberOfEvent();
  if ( nofEvents == 0 ) return;
  
  // print histogram statistics
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->GetH1(1) ) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
       << " Emylar : mean = " << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy") 
               << " rms = " << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") 
               << G4endl;
    G4cout 	       
       << " EGap : mean = " << G4BestUnit(analysisManager->GetH1(2)->mean(), "Energy") 
               << " rms = " << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Energy") 
               << G4endl;
    G4cout
       << " Eelectrode : mean = " << G4BestUnit(analysisManager->GetH1(3)->mean(), "Energy")
               << " rms = " << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Energy")
               << G4endl;
    G4cout
       << " EGap2 : mean = " << G4BestUnit(analysisManager->GetH1(4)->mean(), "Energy")
               << " rms = " << G4BestUnit(analysisManager->GetH1(4)->rms(),  "Energy")
               << G4endl;
    G4cout
       << " Esurface : mean = " << G4BestUnit(analysisManager->GetH1(5)->mean(), "Energy")
               << " rms = " << G4BestUnit(analysisManager->GetH1(5)->rms(),  "Energy")
               << G4endl;
    G4cout
       << " Lmylar : mean = " << G4BestUnit(analysisManager->GetH1(6)->mean(), "Length")
               << " rms = " << G4BestUnit(analysisManager->GetH1(6)->rms(),  "Length")
               << G4endl;
    G4cout 
       << " LGap : mean = " << G4BestUnit(analysisManager->GetH1(7)->mean(), "Length") 
               << " rms = " << G4BestUnit(analysisManager->GetH1(7)->rms(),  "Length") 
               << G4endl;
    G4cout 
       << " Lelectrode : mean = " << G4BestUnit(analysisManager->GetH1(8)->mean(), "Length") 
               << " rms = " << G4BestUnit(analysisManager->GetH1(8)->rms(),  "Length") 
               << G4endl;
    G4cout
       << " LGap2 : mean = " << G4BestUnit(analysisManager->GetH1(9)->mean(), "Length")
               << " rms = " << G4BestUnit(analysisManager->GetH1(9)->rms(),  "Length")
               << G4endl;
    G4cout
       << " Lsurface : mean = " << G4BestUnit(analysisManager->GetH1(10)->mean(), "Length")
               << " rms = " << G4BestUnit(analysisManager->GetH1(10)->rms(),  "Length")
               << G4endl;
  }
  
  // save histograms 
  //
  analysisManager->Write();
  analysisManager->CloseFile();
  
  // complete cleanup
  //
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
