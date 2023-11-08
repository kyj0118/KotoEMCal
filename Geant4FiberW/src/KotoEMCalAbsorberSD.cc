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
/// \file KotoEMCalAbsorberSD.cc
/// \brief Implementation of the KotoEMCalAbsorberSD class

// This project class
#include "KotoEMCalAbsorberHit.hh"
#include "KotoEMCalAbsorberSD.hh"

// Geant4 class
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KotoEMCalAbsorberSD::KotoEMCalAbsorberSD(G4String name, G4int nLayers, G4int nModulesXY)
    : G4VSensitiveDetector(name), fNameSD(name), fNumberOfLayers(nLayers), fNumberOfModulesXY(nModulesXY), fHitsCollection(nullptr), fHCID(-1) {
  fEdep.clear();
  fEweightedx.clear();
  fEweightedy.clear();
  fEweightedz.clear();
  fEweightedt.clear();

  fNumberOfTotalModules = fNumberOfLayers * fNumberOfModulesXY;
  fEdep.resize(fNumberOfTotalModules, 0);
  fEweightedx.resize(fNumberOfTotalModules, 0);
  fEweightedy.resize(fNumberOfTotalModules, 0);
  fEweightedz.resize(fNumberOfTotalModules, 0);
  fEweightedt.resize(fNumberOfTotalModules, 0);

  G4String nameOfHC = name + "HitCollection";
  collectionName.insert(nameOfHC);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KotoEMCalAbsorberSD::~KotoEMCalAbsorberSD() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KotoEMCalAbsorberSD::Initialize(G4HCofThisEvent* hce) {
  fHitsCollection = new KotoEMCalAbsorberHitsCollection(fNameSD, collectionName[0]);
  if (fHCID < 0) {
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }
  hce->AddHitsCollection(fHCID, fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool KotoEMCalAbsorberSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
  auto edep = step->GetTotalEnergyDeposit();
  if (edep != 0) {
    auto prepoint = step->GetPreStepPoint();
    auto postpoint = step->GetPostStepPoint();
    auto touchable = prepoint->GetTouchable();
    G4int copyID = touchable->GetReplicaNumber();
    G4int moduleID = touchable->GetReplicaNumber(1);
    G4double prex = (prepoint->GetPosition()).x();
    G4double prey = (prepoint->GetPosition()).y();
    G4double prez = (prepoint->GetPosition()).z();
    G4double pret = prepoint->GetGlobalTime();

    G4double postx = (postpoint->GetPosition()).x();
    G4double posty = (postpoint->GetPosition()).y();
    G4double postz = (postpoint->GetPosition()).z();
    G4double postt = postpoint->GetGlobalTime();

    G4double x = (prex + postx) / 2.0;
    G4double y = (prey + posty) / 2.0;
    G4double z = (prez + postz) / 2.0;
    G4double t = (pret + postt) / 2.0;

    fEdep[moduleID] += edep;
    fEweightedx[moduleID] += x * edep;
    fEweightedy[moduleID] += y * edep;
    fEweightedz[moduleID] += z * edep;
    fEweightedt[moduleID] += t * edep;
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KotoEMCalAbsorberSD::EndOfEvent(G4HCofThisEvent* hce) {
  for (int i = 0; i < fNumberOfTotalModules; i++) {
    if (fEdep[i] == 0) continue;

    fEweightedx[i] /= fEdep[i];
    fEweightedy[i] /= fEdep[i];
    fEweightedz[i] /= fEdep[i];
    fEweightedt[i] /= fEdep[i];

    fHitsCollection->insert(new KotoEMCalAbsorberHit(fHCID));
    G4int CurrentHitID = fHitsCollection->GetSize() - 1;
    auto hit = (KotoEMCalAbsorberHit*)((hce->GetHC(fHCID))->GetHit(CurrentHitID));

    hit->SetXYZTE(fEweightedx[i], fEweightedy[i], fEweightedz[i], fEweightedt[i], fEdep[i]);
    G4int layerID = i / fNumberOfModulesXY;
    G4int moduleXYID = i % fNumberOfModulesXY;
    hit->SetLayerID(layerID);
    hit->SetSegmentID(moduleXYID);
    hit->SetModuleID(i);
    fEdep[i] = 0.0;
    fEweightedx[i] = 0.0;
    fEweightedy[i] = 0.0;
    fEweightedz[i] = 0.0;
    fEweightedt[i] = 0.0;
  }
}
