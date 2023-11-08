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
/// \file KotoEMCalEventAction.cc
/// \brief Implementation of the KotoEMCalEventAction class

// Root class
#include <Riostream.h>

#include "TClonesArray.h"
#include "TObject.h"
#include "TTree.h"

// This project class
#include "KotoEMCalAbsorberHit.hh"
#include "KotoEMCalCsIHit.hh"
#include "KotoEMCalEventAction.hh"
#include "KotoEMCalPrimaryGeneratorAction.hh"
#include "KotoEMCalRunAction.hh"
#include "KotoEMCalScintillatorHit.hh"
#include "KotoEMCalTriggerCounterHit.hh"

// Geant4 class
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VHitsCollection.hh"
#include "G4ios.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern G4ThreeVector gPrimaryParticlePosition;
extern G4ThreeVector gPrimaryParticleMomentumDirection;
extern int gPrimaryParticlePDG;
extern double gPrimaryParticleEnergy;
extern double gPrimaryParticleMass;
KotoEMCalEventAction::KotoEMCalEventAction(KotoEMCalRunAction *runAction, TTree *tr)
    : G4UserEventAction(), fRunAction(runAction), fTree(tr) {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KotoEMCalEventAction::~KotoEMCalEventAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KotoEMCalEventAction::BeginOfEventAction(const G4Event *) {
  SetRunID(fRunAction->GetRunID());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KotoEMCalEventAction::EndOfEventAction(const G4Event *event) {
  EventInfo.eventID = event->GetEventID();
  if (EventInfo.eventID % 100 == 0) G4cout << "Event No. : " << EventInfo.eventID << "/" << fRunAction->GetNumberOfEvents() << G4endl;

  auto runManager = G4RunManager::GetRunManager();
  auto primaryGenerator = (KotoEMCalPrimaryGeneratorAction *)runManager->GetUserPrimaryGeneratorAction();
  auto particleSource = primaryGenerator->GetParticleSource();

  auto hce = event->GetHCofThisEvent();
  int iarrayScintHit = 0;
  int iarrayX0Hit = 0;
  int iarrayTriggerHit = 0;
  int iarrayAbsorberHit = 0;
  int iarrayCsIHit = 0;
  for (int i = 0; i < hce->GetCapacity(); i++) {
    if (hce->GetHC(i)->GetSize() == 0) continue;
    G4String iHCName = hce->GetHC(i)->GetName();

    // Scintillator Hit
    if (iHCName == "ScintSDHitCollection") {
      int nScintHits = hce->GetHC(i)->GetSize();
      for (int ih = 0; ih < nScintHits; ih++) {
        auto hit = (KotoEMCalScintillatorHit *)(hce->GetHC(i)->GetHit(ih));
        ScintHit.cid[iarrayScintHit] = hit->GetModuleID();
        ScintHit.lid[iarrayScintHit] = hit->GetLayerID();
        ScintHit.segid[iarrayScintHit] = hit->GetSegmentID();
        double xx, yy, zz, tt, ee;
        hit->GetXYZTE(xx, yy, zz, tt, ee);

        ScintHit.x[iarrayScintHit] = xx;
        ScintHit.y[iarrayScintHit] = yy;
        ScintHit.z[iarrayScintHit] = zz;
        ScintHit.t[iarrayScintHit] = tt;
        ScintHit.e[iarrayScintHit] = ee;
        iarrayScintHit++;
      }
    }

    // X0Finder Hit
    else if (iHCName == "X0FinderSDHitCollection") {
      int nX0Hits = hce->GetHC(i)->GetSize();
      for (int ih = 0; ih < nX0Hits; ih++) {
        auto hit = (KotoEMCalTriggerCounterHit *)(hce->GetHC(i)->GetHit(ih));
        X0FinderHit.cid[iarrayX0Hit] = hit->GetModuleID();
        X0FinderHit.lid[iarrayX0Hit] = hit->GetLayerID();
        X0FinderHit.segid[iarrayX0Hit] = hit->GetSegmentID();
        double xx, yy, zz, tt, ee;
        hit->GetXYZTE(xx, yy, zz, tt, ee);

        X0FinderHit.x[iarrayX0Hit] = xx;
        X0FinderHit.y[iarrayX0Hit] = yy;
        X0FinderHit.z[iarrayX0Hit] = zz;
        X0FinderHit.t[iarrayX0Hit] = tt;
        X0FinderHit.e[iarrayX0Hit] = ee;
        iarrayX0Hit++;
      }
    }

    // X0Finder Hit
    else if (iHCName == "TriggerCounterSDHitCollection") {
      int nTriggerHits = hce->GetHC(i)->GetSize();
      for (int ih = 0; ih < nTriggerHits; ih++) {
        auto hit = (KotoEMCalTriggerCounterHit *)(hce->GetHC(i)->GetHit(ih));
        TriggerCounterHit.cid[iarrayTriggerHit] = hit->GetModuleID();
        TriggerCounterHit.lid[iarrayTriggerHit] = hit->GetLayerID();
        TriggerCounterHit.segid[iarrayTriggerHit] = hit->GetSegmentID();
        double xx, yy, zz, tt, ee;
        hit->GetXYZTE(xx, yy, zz, tt, ee);

        TriggerCounterHit.x[iarrayTriggerHit] = xx;
        TriggerCounterHit.y[iarrayTriggerHit] = yy;
        TriggerCounterHit.z[iarrayTriggerHit] = zz;
        TriggerCounterHit.t[iarrayTriggerHit] = tt;
        TriggerCounterHit.e[iarrayTriggerHit] = ee;
        iarrayTriggerHit++;
      }
    }

    // Absorber Hit
    else if (iHCName == "AbsorberSDHitCollection") {
      int nAbsorberHits = hce->GetHC(i)->GetSize();
      for (int ih = 0; ih < nAbsorberHits; ih++) {
        auto hit = (KotoEMCalAbsorberHit *)(hce->GetHC(i)->GetHit(ih));
        AbsorberHit.cid[iarrayAbsorberHit] = hit->GetModuleID();
        AbsorberHit.lid[iarrayAbsorberHit] = hit->GetLayerID();
        AbsorberHit.segid[iarrayAbsorberHit] = hit->GetSegmentID();
        double xx, yy, zz, tt, ee;
        hit->GetXYZTE(xx, yy, zz, tt, ee);
        
        AbsorberHit.x[iarrayAbsorberHit] = xx;
        AbsorberHit.y[iarrayAbsorberHit] = yy;
        AbsorberHit.z[iarrayAbsorberHit] = zz;
        AbsorberHit.t[iarrayAbsorberHit] = tt;
        AbsorberHit.e[iarrayAbsorberHit] = ee;
        iarrayAbsorberHit++;
      }
    }

    // CsI Hit
    else if (iHCName == "CsISDHitCollection") {
      auto hit = (KotoEMCalCsIHit *)(hce->GetHC(i)->GetHit(0));
      CsIHit.e[iarrayCsIHit] = hit->GetEdep();
      CsIHit.xid[iarrayCsIHit] = hit->GetXID();
      CsIHit.yid[iarrayCsIHit] = hit->GetYID();
      CsIHit.cid[iarrayCsIHit] = hit->GetCellID();
      iarrayCsIHit++;
    }
  }

  ScintHit.nhit = iarrayScintHit;
  AbsorberHit.nhit = iarrayAbsorberHit;
  CsIHit.nhit = iarrayCsIHit;
  X0FinderHit.nhit = iarrayX0Hit;
  TriggerCounterHit.nhit = iarrayTriggerHit;

  PrimaryParticle.x = gPrimaryParticlePosition.getX();
  PrimaryParticle.y = gPrimaryParticlePosition.getY();
  PrimaryParticle.z = gPrimaryParticlePosition.getZ();

  PrimaryParticle.e = gPrimaryParticleEnergy;
  PrimaryParticle.m = gPrimaryParticleMass;
  PrimaryParticle.p = sqrt(PrimaryParticle.e * PrimaryParticle.e - PrimaryParticle.m * PrimaryParticle.m);  // magnitude of momentum

  PrimaryParticle.px = gPrimaryParticleMomentumDirection.getX() * PrimaryParticle.p;
  PrimaryParticle.py = gPrimaryParticleMomentumDirection.getY() * PrimaryParticle.p;
  PrimaryParticle.pz = gPrimaryParticleMomentumDirection.getZ() * PrimaryParticle.p;
  fTree->Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KotoEMCalEventAction::SetBranch() {
  fTree->Branch("eventID", &EventInfo.eventID, "eventID/I");
  fTree->Branch("runID", &EventInfo.runID, "runID/I");
  fTree->Branch("randomSeed", &EventInfo.randomSeed, "randomSeed/L");

  fTree->Branch("PrimaryParticle.x", &PrimaryParticle.x, "PrimaryParticle.x/D");
  fTree->Branch("PrimaryParticle.y", &PrimaryParticle.y, "PrimaryParticle.y/D");
  fTree->Branch("PrimaryParticle.z", &PrimaryParticle.z, "PrimaryParticle.z/D");
  fTree->Branch("PrimaryParticle.px", &PrimaryParticle.px, "PrimaryParticle.px/D");
  fTree->Branch("PrimaryParticle.py", &PrimaryParticle.py, "PrimaryParticle.py/D");
  fTree->Branch("PrimaryParticle.pz", &PrimaryParticle.pz, "PrimaryParticle.pz/D");
  fTree->Branch("PrimaryParticle.p", &PrimaryParticle.p, "PrimaryParticle.p/D");
  fTree->Branch("PrimaryParticle.m", &PrimaryParticle.m, "PrimaryParticle.m/D");
  fTree->Branch("PrimaryParticle.e", &PrimaryParticle.e, "PrimaryParticle.e/D");
  fTree->Branch("PrimaryParticle.PDG", &PrimaryParticle.PDG, "PrimaryParticle.PDG/I");

  fTree->Branch("nScintHit", &ScintHit.nhit, "nScintHit/I");
  fTree->Branch("ScintHit.ModuleID", ScintHit.cid, "ScintHit.ModuleID[nScintHit]/I");
  fTree->Branch("ScintHit.LayerID", ScintHit.lid, "ScintHit.LayerID[nScintHit]/I");
  fTree->Branch("ScintHit.SegmentID", ScintHit.segid, "ScintHit.SegmentID[nScintHit]/I");
  fTree->Branch("ScintHit.x", ScintHit.x, "ScintHit.x[nScintHit]/D");
  fTree->Branch("ScintHit.y", ScintHit.y, "ScintHit.y[nScintHit]/D");
  fTree->Branch("ScintHit.z", ScintHit.z, "ScintHit.z[nScintHit]/D");
  fTree->Branch("ScintHit.t", ScintHit.t, "ScintHit.t[nScintHit]/D");
  fTree->Branch("ScintHit.e", ScintHit.e, "ScintHit.e[nScintHit]/D");

  fTree->Branch("nX0FinderHit", &X0FinderHit.nhit, "nX0FinderHit/I");
  fTree->Branch("X0FinderHit.ModuleID", X0FinderHit.cid, "X0FinderHit.ModuleID[nX0FinderHit]/I");
  fTree->Branch("X0FinderHit.LayerID", X0FinderHit.lid, "X0FinderHit.LayerID[nX0FinderHit]/I");
  fTree->Branch("X0FinderHit.SegmentID", X0FinderHit.segid, "X0FinderHit.SegmentID[nX0FinderHit]/I");
  fTree->Branch("X0FinderHit.x", X0FinderHit.x, "X0FinderHit.x[nX0FinderHit]/D");
  fTree->Branch("X0FinderHit.y", X0FinderHit.y, "X0FinderHit.y[nX0FinderHit]/D");
  fTree->Branch("X0FinderHit.z", X0FinderHit.z, "X0FinderHit.z[nX0FinderHit]/D");
  fTree->Branch("X0FinderHit.t", X0FinderHit.t, "X0FinderHit.t[nX0FinderHit]/D");
  fTree->Branch("X0FinderHit.e", X0FinderHit.e, "X0FinderHit.e[nX0FinderHit]/D");

  fTree->Branch("nTriggerCounterHit", &TriggerCounterHit.nhit, "nTriggerCounterHit/I");
  fTree->Branch("TriggerCounterHit.ModuleID", TriggerCounterHit.cid, "TriggerCounterHit.ModuleID[nTriggerCounterHit]/I");
  fTree->Branch("TriggerCounterHit.LayerID", TriggerCounterHit.lid, "TriggerCounterHit.LayerID[nTriggerCounterHit]/I");
  fTree->Branch("TriggerCounterHit.SegmentID", TriggerCounterHit.segid, "TriggerCounterHit.SegmentID[nTriggerCounterHit]/I");
  fTree->Branch("TriggerCounterHit.x", TriggerCounterHit.x, "TriggerCounterHit.x[nTriggerCounterHit]/D");
  fTree->Branch("TriggerCounterHit.y", TriggerCounterHit.y, "TriggerCounterHit.y[nTriggerCounterHit]/D");
  fTree->Branch("TriggerCounterHit.z", TriggerCounterHit.z, "TriggerCounterHit.z[nTriggerCounterHit]/D");
  fTree->Branch("TriggerCounterHit.t", TriggerCounterHit.t, "TriggerCounterHit.t[nTriggerCounterHit]/D");
  fTree->Branch("TriggerCounterHit.e", TriggerCounterHit.e, "TriggerCounterHit.e[nTriggerCounterHit]/D");

  fTree->Branch("nAbsorberHit", &AbsorberHit.nhit, "nAbsorberHit/I");
  fTree->Branch("AbsorberHit.ModuleID", AbsorberHit.cid, "AbsorberHit.ModuleID[nAbsorberHit]/I");
  fTree->Branch("AbsorberHit.LayerID", AbsorberHit.lid, "AbsorberHit.LayerID[nAbsorberHit]/I");
  fTree->Branch("AbsorberHit.SegmentID", AbsorberHit.segid, "AbsorberHit.SegmentID[nAbsorberHit]/I");
  fTree->Branch("AbsorberHit.x", AbsorberHit.x, "AbsorberHit.x[nAbsorberHit]/D");
  fTree->Branch("AbsorberHit.y", AbsorberHit.y, "AbsorberHit.y[nAbsorberHit]/D");
  fTree->Branch("AbsorberHit.z", AbsorberHit.z, "AbsorberHit.z[nAbsorberHit]/D");
  fTree->Branch("AbsorberHit.t", AbsorberHit.t, "AbsorberHit.t[nAbsorberHit]/D");
  fTree->Branch("AbsorberHit.e", AbsorberHit.e, "AbsorberHit.e[nAbsorberHit]/D");

  fTree->Branch("nCsIHit", &CsIHit.nhit, "nCsIHit/I");
  fTree->Branch("CsIHit.CellID", CsIHit.cid, "CsIHit.CellID[nCsIHit]/I");
  fTree->Branch("CsIHit.xid", CsIHit.xid, "CsIHit.xid[nCsIHit]/I");
  fTree->Branch("CsIHit.yid", CsIHit.yid, "CsIHit.yid[nCsIHit]/I");
  fTree->Branch("CsIHit.e", CsIHit.e, "CsIHit.e[nCsIHit]/D");
}
void KotoEMCalEventAction::SetRunID(G4int RunID) {
  EventInfo.runID = RunID;
}
void KotoEMCalEventAction::SetRandomSeed(long seed) {
  EventInfo.randomSeed = seed;
}
