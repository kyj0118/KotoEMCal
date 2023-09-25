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
/// \file KotoEMCalCsISD.hh
/// \brief Definition of the KotoEMCalCsISD class

#ifndef KotoEMCalCsISD_h
#define KotoEMCalCsISD_h 1

#include "G4VSensitiveDetector.hh"
#include "KotoEMCalCsIHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

/// EM calorimeter sensitive detector
using namespace std;
class KotoEMCalCsISD : public G4VSensitiveDetector
{   
public:
  KotoEMCalCsISD(G4String name, G4int xid, G4int yid);
  virtual ~KotoEMCalCsISD();
  
  virtual void Initialize(G4HCofThisEvent*HCE);
  virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
  void EndOfEvent(G4HCofThisEvent*);
private:
  G4String fNameSD;
  G4int fLayerId;
  G4int fSegmentId;

  G4int fxid;
  G4int fyid;

  G4double fEdep;
  
  KotoEMCalCsIHitsCollection* fHitsCollection;
  G4int fHCID;
  
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif