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
/// \file PPACEventAction.cc
/// \brief Implementation of the PPACEventAction class

#include "PPACEventAction.hh"
#include "PPACEventActionMessenger.hh"
#include "PPACCalorimeterSD.hh"
#include "PPACCalorHit.hh"
#include "B4Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PPACEventAction::PPACEventAction()
 : G4UserEventAction(),
   fMessenger(0),
   fPrintModulo(1)
{
  fMessenger = new PPACEventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PPACEventAction::~PPACEventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PPACCalorHitsCollection* 
PPACEventAction::GetHitsCollection(const G4String& hcName,
                                  const G4Event* event) const
{
  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(hcName);
  PPACCalorHitsCollection* hitsCollection 
    = static_cast<PPACCalorHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4cerr << "Cannot access hitsCollection " << hcName << G4endl;
    exit(1);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PPACEventAction::PrintEventStatistics(
                              G4double mylarEdep, G4double mylarTrackLength,
                              G4double gapEdep, G4double gapTrackLength,
			      G4double electroEdep, G4double electroTrackLength,
			      G4double gap2Edep, G4double gap2TrackLength,
			      G4double surfaceEdep, G4double surfaceTrackLength) const
{
  // print event statistics
  G4cout
     << "   mylar: total energy: " 
     << std::setw(7) << G4BestUnit(mylarEdep, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(mylarTrackLength, "Length")
     << G4endl
     << "        Gap: total energy: " 
     << std::setw(7) << G4BestUnit(gapEdep, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(gapTrackLength, "Length")
     << G4endl
     << "   Electrode: total energy: "
     << std::setw(7) << G4BestUnit(electroEdep, "Energy")
     << "       total track length: "
     << std::setw(7) << G4BestUnit(electroTrackLength, "Length")
     << G4endl
     << "        Gap2: total energy: "
     << std::setw(7) << G4BestUnit(gap2Edep, "Energy")
     << "       total track length: "
     << std::setw(7) << G4BestUnit(gap2TrackLength, "Length")
     << G4endl
     << "       Surface: total energy: "
     << std::setw(7) << G4BestUnit(surfaceEdep, "Energy")
     << "       total track length: "
     << std::setw(7) << G4BestUnit(surfaceTrackLength, "Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PPACEventAction::BeginOfEventAction(const G4Event* event)
{  

  G4int eventID = event->GetEventID();
  if ( eventID % fPrintModulo == 0 )  { 
    G4cout << "\n---> Begin of event: " << eventID << G4endl;
    //CLHEP::HepRandom::showEngineStatus();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PPACEventAction::EndOfEventAction(const G4Event* event)
{  
//  G4HCofThisEvent* HCTE = event-> GetHCofThisEvent();
//  if (!HCTE) return;  // no hits in this events. nothing to do!

//  G4THitsCollection<MylarHitsCollection> *CHC_TPC = NULL;


  // Get hits collections
  PPACCalorHitsCollection* mylarHC
    = GetHitsCollection("MylarHitsCollection", event);
  PPACCalorHitsCollection* gapHC
    = GetHitsCollection("GapHitsCollection", event);
  PPACCalorHitsCollection* electroHC
    = GetHitsCollection("ElectroHitsCollection", event);
  PPACCalorHitsCollection* gap2HC
    = GetHitsCollection("Gap2HitsCollection", event);
  PPACCalorHitsCollection* surfaceHC
    = GetHitsCollection("SurfaceHitsCollection", event);
  // Get hit with total values
  PPACCalorHit* mylarHit = (*mylarHC)[mylarHC->entries()-1];
  PPACCalorHit* gapHit = (*gapHC)[gapHC->entries()-1];
  PPACCalorHit* electroHit = (*electroHC)[electroHC->entries()-1];
  PPACCalorHit* gap2Hit = (*gap2HC)[gap2HC->entries()-1];
  PPACCalorHit* surfaceHit = (*surfaceHC)[surfaceHC->entries()-1];
  // Print per event (modulo n)
  //
  G4int eventID = event->GetEventID();
  if ( eventID % fPrintModulo == 0) {
    G4cout << "---> End of event: " << eventID << G4endl;     

    PrintEventStatistics(
      mylarHit->GetEdep(), mylarHit->GetTrackLength(),
      gapHit->GetEdep(), gapHit->GetTrackLength(),
      electroHit->GetEdep(), electroHit->GetTrackLength(),
      gap2Hit->GetEdep(), gap2Hit->GetTrackLength(),
      surfaceHit->GetEdep(), surfaceHit->GetTrackLength());
  }  
  
  // Fill histograms, ntuple
  //
//  G4int nHits = mylarHC -> entries();

//  G4ThreeVector prePos = (mylarHit)[1]->GetPrePosition();
  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 
  // fill histograms
  analysisManager->FillH1(1, mylarHit->GetEdep());
  analysisManager->FillH1(2, gapHit->GetEdep());
  analysisManager->FillH1(3, electroHit->GetEdep());
  analysisManager->FillH1(4, gap2Hit->GetEdep());
  analysisManager->FillH1(5, surfaceHit->GetEdep());
  analysisManager->FillH1(6, mylarHit->GetTrackLength());
  analysisManager->FillH1(7, gapHit->GetTrackLength());
  analysisManager->FillH1(8, electroHit->GetTrackLength());
  analysisManager->FillH1(9, gap2Hit->GetTrackLength());
  analysisManager->FillH1(10, surfaceHit->GetTrackLength());
  // fill ntuple
  analysisManager->FillNtupleDColumn(0, mylarHit->GetEdep());
  analysisManager->FillNtupleDColumn(1, gapHit->GetEdep());
  analysisManager->FillNtupleDColumn(2, electroHit->GetEdep());
  analysisManager->FillNtupleDColumn(3, gap2Hit->GetEdep());
  analysisManager->FillNtupleDColumn(4, surfaceHit->GetEdep());
  analysisManager->FillNtupleDColumn(5, gap2Hit->GetPrePosition().x());
  analysisManager->FillNtupleDColumn(6, gap2Hit->GetPostPosition().x());
  analysisManager->FillNtupleDColumn(7, gap2Hit->GetPrePosition().z());
  analysisManager->FillNtupleDColumn(8, gap2Hit->GetPostPosition().z());
  analysisManager->FillNtupleDColumn(9, surfaceHit->GetTrackLength());
  analysisManager->AddNtupleRow();  

/*  for (G4int i = 0; i != nHits; i++) {
   G4ThreeVector prePos = (*mylarHC)[i]->GetPrePosition();
   G4ThreeVector postPos = (*mylarHC)[i]->GetPostPosition();
   analysisManager->FillNtupleDColumn(0, 0.);
   analysisManager->FillNtupleDColumn(5, prePos.z());
   analysisManager->FillNtupleDColumn(6, postPos.z());
   G4cout << "prePos.z = " << prePos.z() <<" "<<prePos.x()<< G4endl;
   G4cout << "postPos.z = " << postPos.z() << " "<<postPos.x() <<G4endl;
//   G4ThreeVector postPos2 = (*gap2HC)[i]->GetPostPosition();
//   G4cout << "postPos2.x = " << postPos2.x() << G4endl;
  }*/
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
