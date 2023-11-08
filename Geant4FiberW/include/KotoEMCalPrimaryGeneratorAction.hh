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
/// \file KotoEMCalPrimaryGeneratorAction.hh
/// \brief Definition of the KotoEMCalPrimaryGeneratorAction class

#ifndef KotoEMCalPrimaryGeneratorAction_h
#define KotoEMCalPrimaryGeneratorAction_h 1

// Geant4 class
#include "G4GeneralParticleSource.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4GenericMessenger;
class G4Event;
class G4ParticleDefinition;

class KotoEMCalPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  KotoEMCalPrimaryGeneratorAction();
  virtual ~KotoEMCalPrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event*);

  void SetOptionConfiguration(std::map<std::string, std::string> optConfig);
  G4VPrimaryGenerator* GetParticleSource();

  std::map<std::string, std::string> fOptionConfiguration;
  G4String fEnergyOrMomentum;
  G4double fEnergyMin;
  G4double fEnergyMax;
  G4int fEnergyNstep;
  G4double fEnergyStepSize;

  G4double fTheta;
  G4double fThetaMin;
  G4double fThetaMax;
  G4int fThetaNstep;
  G4double fThetaStepSize;

  G4double fPhi;
  G4double fPhiMin;
  G4double fPhiMax;
  G4int fPhiNstep;
  G4double fPhiStepSize;

  G4double fDirX;
  G4double fDirY;
  G4double fDirZ;

  G4double fPosX;
  G4double fPosY;
  G4double fPosZ;
  G4bool fPosRandomX;
  G4bool fPosRandomY;
  G4bool fPosRandomZ;
  G4double fMinX;
  G4double fMinY;
  G4double fMinZ;
  G4double fMaxX;
  G4double fMaxY;
  G4double fMaxZ;

  G4String fGenerateEnergyOption;
  G4String fGenerateAngleOption;
  G4String fGeneratePositionOption;

  G4String fParticleName;
  G4ParticleGun* fParticleGun;
  G4GeneralParticleSource* fGeneralParticleSource;
  G4ParticleDefinition* fParticle;
};

#endif
