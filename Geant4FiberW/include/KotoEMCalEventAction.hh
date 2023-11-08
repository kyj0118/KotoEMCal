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
/// \file KotoEMCalEventAction.hh
/// \brief Definition of the KotoEMCalEventAction class

#ifndef KotoEMCalEventAction_h
#define KotoEMCalEventAction_h 1

// Root class
#include "TInterpreter.h"
#include "TSystem.h"
#include "TTree.h"

// c++ std
#include <array>
#include <vector>

// This project class
#include "KotoEMCalRunAction.hh"

// Geant4 class
#include "G4UserEventAction.hh"
#include "globals.hh"

using namespace std;
/// Event action
const int kMaxLayer = 24;
const int kMaxSegment = 16;
const int kMaxScintillator = kMaxLayer * kMaxSegment;

struct ScintHitStruct {
  int nhit;
  int nMaximumHits = kMaxScintillator;
  int cid[kMaxScintillator];
  int lid[kMaxScintillator];
  int segid[kMaxScintillator];
  double x[kMaxScintillator];
  double y[kMaxScintillator];
  double z[kMaxScintillator];
  double t[kMaxScintillator];
  double e[kMaxScintillator];
};

struct AbsorberHitStruct {
  int nhit;
  int nMaximumHits = kMaxScintillator;
  int cid[kMaxScintillator];
  int lid[kMaxScintillator];
  int segid[kMaxScintillator];
  double x[kMaxScintillator];
  double y[kMaxScintillator];
  double z[kMaxScintillator];
  double t[kMaxScintillator];
  double e[kMaxScintillator];
};

struct CsIHitStruct {
  int nhit;
  int cid[99];
  int xid[99];
  int yid[99];
  double e[99];
};

struct EventInfoStruct {
  int eventID;
  int runID;
  long randomSeed;
};

struct PrimaryParticleInfoStruct {
  double x, y, z;
  double px, py, pz;
  double p, m, e;
  int PDG;
};

class KotoEMCalEventAction : public G4UserEventAction {
 public:
  KotoEMCalEventAction(KotoEMCalRunAction *runAction, TTree *tr);
  virtual ~KotoEMCalEventAction();
  void SetBranch();
  void SetRunID(G4int RunID);
  void SetRandomSeed(long seed);
  virtual void BeginOfEventAction(const G4Event *);
  virtual void EndOfEventAction(const G4Event *);

  ScintHitStruct ScintHit;
  ScintHitStruct X0FinderHit;
  ScintHitStruct TriggerCounterHit;
  AbsorberHitStruct AbsorberHit;
  CsIHitStruct CsIHit;
  EventInfoStruct EventInfo;
  PrimaryParticleInfoStruct PrimaryParticle;

 private:
  KotoEMCalRunAction *fRunAction;
  TTree *fTree;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
