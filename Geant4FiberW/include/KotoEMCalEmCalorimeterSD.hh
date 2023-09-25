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
/// \file KotoEMCalEmCalorimeterSD.hh
/// \brief Definition of the KotoEMCalEmCalorimeterSD class

#ifndef KotoEMCalEmCalorimeterSD_h
#define KotoEMCalEmCalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "KotoEMCalEmCalorimeterHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

/// EM calorimeter sensitive detector
using namespace std;
class KotoEMCalEmCalorimeterSD : public G4VSensitiveDetector
{   
public:
  KotoEMCalEmCalorimeterSD(G4String name, G4int layerNumber, G4int scintNumber);
  virtual ~KotoEMCalEmCalorimeterSD();
  
  virtual void Initialize(G4HCofThisEvent*HCE);
  virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
  void EndOfEvent(G4HCofThisEvent*);
private:
  G4String fNameSD;
  G4int fLayerId;
  G4int fSegmentId;

  vector<G4double> fEdep;
  vector<G4double> fEweightedx;
  vector<G4double> fEweightedy;
  vector<G4double> fEweightedz;
  vector<G4double> fEweightedt;

  vector<vector<G4double>> fStepEdep;
  
  vector<vector<G4double>> fPreStepx;
  vector<vector<G4double>> fPreStepy;
  vector<vector<G4double>> fPreStepz;
  vector<vector<G4double>> fPreStept;
  
  vector<vector<G4double>> fPostStepx;
  vector<vector<G4double>> fPostStepy;
  vector<vector<G4double>> fPostStepz;
  vector<vector<G4double>> fPostStept;
  
  vector<vector<G4double>> fParticlePx;
  vector<vector<G4double>> fParticlePy;
  vector<vector<G4double>> fParticlePz;
  vector<vector<G4int>> fParticleTrackID;
  vector<vector<G4int>> fParticleParentID;
  vector<vector<G4double>> fParticleCharge;
  vector<vector<G4double>> fParticleMass;
  vector<vector<G4int>> fParticlePDGID;
  
  KotoEMCalEmCalorimeterHitsCollection* fHitsCollection;
  G4int fHCID;
  
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
