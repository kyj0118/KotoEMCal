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

#ifndef KotoEMCalLeadSD_h
#define KotoEMCalLeadSD_h 1

// This project class
#include "KotoEMCalLeadHit.hh"

// Geant4 class
#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

/// EM calorimeter sensitive detector

class KotoEMCalLeadSD : public G4VSensitiveDetector {
 public:
  KotoEMCalLeadSD(G4String name, G4int layerNumber);
  virtual ~KotoEMCalLeadSD();

  virtual void Initialize(G4HCofThisEvent* HCE);
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  void EndOfEvent(G4HCofThisEvent*);

 private:
  G4String fNameSD;
  G4int fLayerId;
  G4double fEdep;
  G4double fEweightedx;
  G4double fEweightedy;
  G4double fEweightedz;
  G4double fEweightedt;
  std::vector<G4double> fStepEdep;

  std::vector<G4double> fPreStepx;
  std::vector<G4double> fPreStepy;
  std::vector<G4double> fPreStepz;
  std::vector<G4double> fPreStept;

  std::vector<G4double> fPostStepx;
  std::vector<G4double> fPostStepy;
  std::vector<G4double> fPostStepz;
  std::vector<G4double> fPostStept;

  std::vector<G4double> fParticlePx;
  std::vector<G4double> fParticlePy;
  std::vector<G4double> fParticlePz;
  std::vector<G4int> fParticleTrackID;
  std::vector<G4int> fParticleParentID;
  std::vector<G4double> fParticleCharge;
  std::vector<G4double> fParticleMass;
  std::vector<G4int> fParticlePDGID;

  KotoEMCalLeadHitsCollection* fHitsCollection;
  G4int fHCID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
