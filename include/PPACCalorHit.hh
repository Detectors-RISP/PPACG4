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
/// \file PPACCalorHit.hh
/// \brief Definition of the PPACCalorHit class

#ifndef PPACCalorHit_h
#define PPACCalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

/// Calorimeter hit class
///
/// It defines data members to store the the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

class PPACCalorHit : public G4VHit
{
  public:
    PPACCalorHit();
    PPACCalorHit(const PPACCalorHit&);
    virtual ~PPACCalorHit();

    // operators
    const PPACCalorHit& operator=(const PPACCalorHit&);
    G4int operator==(const PPACCalorHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    void Add(G4double de, G4double dl, G4ThreeVector dpre, G4ThreeVector dpost);

    // get methods
    G4double GetEdep() const;
    G4double GetTrackLength() const;
    const G4ThreeVector &GetPrePosition() const;     
    const G4ThreeVector &GetPostPosition() const;
  
  private:
    G4double fEdep;        ///< Energy deposit in the sensitive volume
    G4double fTrackLength; ///< Track length in the  sensitive volume
    G4ThreeVector fprexyz;
    G4ThreeVector fpostxyz;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<PPACCalorHit> PPACCalorHitsCollection;

extern G4Allocator<PPACCalorHit> PPACCalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* PPACCalorHit::operator new(size_t)
{
  void *hit;
  hit = (void *) PPACCalorHitAllocator.MallocSingle();
  return hit;
}

inline void PPACCalorHit::operator delete(void *hit)
{
  PPACCalorHitAllocator.FreeSingle((PPACCalorHit*) hit);
}

inline void PPACCalorHit::Add(G4double de, G4double dl, G4ThreeVector dpre, G4ThreeVector dpost) {
  fEdep += de; 
  fTrackLength += dl;
  fprexyz += dpre;
  fpostxyz += dpost;
}

inline G4double PPACCalorHit::GetEdep() const { 
  return fEdep; 
}

inline G4double PPACCalorHit::GetTrackLength() const { 
  return fTrackLength; 
}

inline const G4ThreeVector &PPACCalorHit::GetPrePosition() const {
  return fprexyz;
}

inline const G4ThreeVector &PPACCalorHit::GetPostPosition() const {
  return fpostxyz;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
