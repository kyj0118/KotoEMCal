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
/// \file KotoEMCalActionInitialization.cc
/// \brief Implementation of the KotoEMCalActionInitialization class

// This project class
#include "KotoEMCalActionInitialization.hh"
#include "KotoEMCalEventAction.hh"
#include "KotoEMCalPrimaryGeneratorAction.hh"
#include "KotoEMCalRunAction.hh"
#include "KotoEMCalStackingAction.hh"

// Geant4 class
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern bool gSaveStepLevel;

KotoEMCalActionInitialization::KotoEMCalActionInitialization(TTree* tr)
    : G4VUserActionInitialization(), fTree(tr) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KotoEMCalActionInitialization::~KotoEMCalActionInitialization() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KotoEMCalActionInitialization::BuildForMaster() const {
  SetUserAction(new KotoEMCalRunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KotoEMCalActionInitialization::Build() const {
  SetUserAction(new KotoEMCalPrimaryGeneratorAction);
  auto runAction = new KotoEMCalRunAction();
  auto eventAction = new KotoEMCalEventAction(runAction, fTree);
  SetUserAction(eventAction);
  SetUserAction(runAction);

  eventAction->SetRandomSeed(CLHEP::HepRandom::getTheSeed());
  eventAction->SetSaveStepLevel(gSaveStepLevel);
  eventAction->SetBranch();
  // SetUserAction(new KotoEMCalStackingAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
