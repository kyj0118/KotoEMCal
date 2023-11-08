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
/// \file KotoEMCalPrimaryGeneratorAction.cc
/// \brief Implementation of the KotoEMCalPrimaryGeneratorAction class

// Root class
#include "TRandom3.h"

// This project class
#include "KotoEMCalPrimaryGeneratorAction.hh"

// Geant4 class
#include "G4Event.hh"
#include "G4GenericMessenger.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector gPrimaryParticlePosition;
G4ThreeVector gPrimaryParticleMomentumDirection;
int gPrimaryParticlePDG;
double gPrimaryParticleEnergy;
double gPrimaryParticleMass;

using namespace std;

KotoEMCalPrimaryGeneratorAction::KotoEMCalPrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(nullptr),
      fGeneralParticleSource(nullptr),
      fParticle(nullptr) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KotoEMCalPrimaryGeneratorAction::~KotoEMCalPrimaryGeneratorAction() {
  delete fParticleGun;
  delete fGeneralParticleSource;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KotoEMCalPrimaryGeneratorAction::SetOptionConfiguration(map<string, string> config) {
  fOptionConfiguration = config;
  G4bool invalid_flag = false;
  if (config["useGPS"] == "true") {
    fGeneralParticleSource = new G4GeneralParticleSource();
  } else if (config["useGPS"] == "false") {
    auto particleTable = G4ParticleTable::GetParticleTable();
    fParticle = particleTable->FindParticle(config["particle"]);
    if (fParticle) {
      fParticleGun = new G4ParticleGun(fParticle);
    } else {
      cout << "Wrong particle name: " << config["particle"] << endl;
      invalid_flag = true;
    }

    // Energy generation option
    fGenerateEnergyOption = config["energyOpt"];
    if (config["energyOpt"] == "fixed") {
      if (config.find("energy") != config.end()) {
        G4double Ek = atof(config["energy"].data()) * MeV;
        fParticleGun->SetParticleEnergy(Ek);
        fEnergyOrMomentum = "energy";
      } else if (config.find("momentum") != config.end()) {
        G4double mom = atof(config["momentum"].data()) * MeV;
        fParticleGun->SetParticleMomentum(mom);
        fEnergyOrMomentum = "momentum";
      } else {
        cout << "energyOpt=fixed: The option energy or momentum is necessary (in MeV)" << endl;
        invalid_flag = true;
      }
    } else if (config["energyOpt"] == "uniform") {
      if (config.find("energyMin") != config.end() && config.find("energyMax") != config.end()) {
        fEnergyMin = atof(config["energyMin"].data()) * MeV;
        fEnergyMax = atof(config["energyMax"].data()) * MeV;
        fEnergyOrMomentum = "energy";
      } else if (config.find("momentumMin") != config.end() && config.find("momentumMax") != config.end()) {
        fEnergyMin = atof(config["momentumMin"].data()) * MeV;
        fEnergyMax = atof(config["momentumMax"].data()) * MeV;
        fEnergyOrMomentum = "momentum";
      } else {
        cout << "energyOpt=uniform: The options energyMin and energyMax" << endl;
        cout << "(or momentumMin and momentumMax) are necessary (in MeV)" << endl;
        invalid_flag = true;
      }
    } else if (config["energyOpt"] == "step") {
      if (config.find("energyBegin") != config.end() && config.find("energyNstep") != config.end() && config.find("energyStepSize") != config.end()) {
        fEnergyMin = atof(config["energyBegin"].data()) * MeV;
        fEnergyNstep = atoi(config["energyNstep"].data());
        fEnergyStepSize = atof(config["energyStepSize"].data()) * MeV;
        fEnergyOrMomentum = "energy";
      } else if (config.find("momentumBegin") != config.end() && config.find("momentumNstep") != config.end() && config.find("momentumStepSize") != config.end()) {
        fEnergyMin = atof(config["momentumBegin"].data()) * MeV;
        fEnergyNstep = atoi(config["momentumNstep"].data());
        fEnergyStepSize = atof(config["momentumStepSize"].data()) * MeV;
        fEnergyOrMomentum = "momentum";
      } else {
        cout << "energyOpt=step: The options energyBegin, energyNstep and energyStepSize" << endl;
        cout << "(or momentumBegin, momentumNstep and momentumStepSize) are necessary (in MeV)" << endl;
        invalid_flag = true;
      }
    } else {
      cout << "Option energyOpt is necessary [energyOpt=fixed or uniform or step]" << endl;
      invalid_flag = true;
    }

    // Angle generation option
    fGenerateAngleOption = config["angleOpt"];
    if (config["angleOpt"] == "fixed") {
      if (config.find("theta") != config.end() && config.find("phi") != config.end()) {
        fTheta = atof(config["theta"].data()) * deg;
        fPhi = atof(config["phi"].data()) * deg;
        fDirX = sin(fTheta) * cos(fPhi);
        fDirY = sin(fTheta) * sin(fPhi);
        fDirZ = cos(fTheta);
      } else if (config.find("directionX") != config.end() && config.find("directionY") != config.end() && config.find("directionZ") != config.end()) {
        fDirX = atof(config["directionX"].data());
        fDirY = atof(config["directionY"].data());
        fDirZ = atof(config["directionZ"].data());
      } else {
        cout << "angleOpt=fixed: The options theta and phi (or directionX, directionY, directionZ) are necessary" << endl;
        invalid_flag = true;
      }
    } else if (config["angleOpt"] == "uniformSolid") {
      if (config.find("thetaMin") != config.end() && config.find("thetaMax") != config.end()) {
        fThetaMin = atof(config["thetaMin"].data()) * deg;
        fThetaMax = atof(config["thetaMax"].data()) * deg;
        fPhiMin = 0;
        fPhiMax = CLHEP::twopi;
      } else {
        cout << "angleOpt=uniformSolid: The options thetaMin and thetaMax are necessary" << endl;
        invalid_flag = true;
      }
    } else if (config["angleOpt"] == "uniformTheta") {
      if (config.find("thetaMin") != config.end() && config.find("thetaMax") != config.end()) {
        fThetaMin = atof(config["thetaMin"].data()) * deg;
        fThetaMax = atof(config["thetaMax"].data()) * deg;
        fPhiMin = 0;
        fPhiMax = CLHEP::twopi;
      } else {
        cout << "angleOpt=uniformTheta: The options thetaMin and thetaMax are necessary" << endl;
        invalid_flag = true;
      }
    } else if (config["angleOpt"] == "stepTheta") {
      if (config.find("thetaBegin") != config.end() && config.find("thetaNstep") != config.end() && config.find("thetaStepSize") != config.end()) {
        fThetaMin = atof(config["thetaBegin"].data()) * MeV;
        fThetaNstep = atoi(config["thetaNstep"].data());
        fThetaStepSize = atof(config["thetaStepSize"].data()) * MeV;
        fPhiMin = 0;
        fPhiMax = CLHEP::twopi;
      } else {
        cout << "angleOpt=stepTheta: The options thetaBegin, thetaNstep and thetaStepSize are necessary" << endl;
        invalid_flag = true;
      }
    }

    else {
      cout << "Option angleOpt is necessary [angleOpt=fixed or uniformSolid or uniformTheta or stepTheta]" << endl;
      invalid_flag = true;
    }

    // Position generation option
    if (config.find("posX") != config.end()) {
      fPosX = atof(config["posX"].data()) * mm;
      fPosRandomX = false;
    } else if (config.find("minX") != config.end() && config.find("maxX") != config.end()) {
      fMinX = atof(config["minX"].data()) * mm;
      fMaxY = atof(config["maxX"].data()) * mm;
      fPosRandomX = true;
    } else {
      cout << "Primary generation positon X should be defined. Use posX (or minX and maxX) in unit of mm" << endl;
      invalid_flag = true;
    }
    if (config.find("posY") != config.end()) {
      fPosY = atof(config["posY"].data()) * mm;
      fPosRandomY = false;
    } else if (config.find("minY") != config.end() && config.find("maxY") != config.end()) {
      fMinY = atof(config["minY"].data()) * mm;
      fMaxY = atof(config["maxY"].data()) * mm;
      fPosRandomY = true;
    } else {
      cout << "Primary generation positon Y should be defined. Use posY (or minY and maxY) in unit of mm" << endl;
      invalid_flag = true;
    }
    if (config.find("posZ") != config.end()) {
      fPosZ = atof(config["posZ"].data()) * mm;
      fPosRandomZ = false;
    } else if (config.find("minZ") != config.end() && config.find("maxZ") != config.end()) {
      fMinZ = atof(config["minZ"].data()) * mm;
      fMaxZ = atof(config["maxZ"].data()) * mm;
      fPosRandomZ = true;
    } else {
      cout << "Primary generation positon Z should be defined. Use posZ (or minZ and maxZ) in unit of mm" << endl;
      invalid_flag = true;
    }

  } else {
    cout << "Option useGPS is necessary [useGPS=true or false]" << endl;
    invalid_flag = true;
  }

  if (invalid_flag) {
    cerr << "You used an invalid option configuration. Exit" << endl;
    exit(1);
  }
}

G4VPrimaryGenerator* KotoEMCalPrimaryGeneratorAction::GetParticleSource() {
  if (fOptionConfiguration["useGPS"] == "true")
    return fGeneralParticleSource;
  else
    return fParticleGun;
}

void KotoEMCalPrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  if (fOptionConfiguration["useGPS"] == "true") {
    fGeneralParticleSource->GeneratePrimaryVertex(event);
    gPrimaryParticlePosition = fGeneralParticleSource->GetParticlePosition();
    gPrimaryParticleEnergy = fGeneralParticleSource->GetParticleEnergy();
    gPrimaryParticleMomentumDirection = fGeneralParticleSource->GetParticleMomentumDirection();
    gPrimaryParticlePDG = fGeneralParticleSource->GetParticleDefinition()->GetPDGEncoding();
    gPrimaryParticleMass = fGeneralParticleSource->GetParticleDefinition()->GetPDGMass();
  } else {
    fParticleGun->SetParticleDefinition(fParticle);

    // position
    if (fPosRandomX == true) fPosX = gRandom->Uniform(fMinX, fMaxX);
    if (fPosRandomY == true) fPosY = gRandom->Uniform(fMinY, fMaxY);
    if (fPosRandomZ == true) fPosZ = gRandom->Uniform(fMinZ, fMaxZ);
    fParticleGun->SetParticlePosition(G4ThreeVector(fPosX, fPosY, fPosZ));

    // time
    fParticleGun->SetParticleTime(0.0 * ns);

    // energy
    if (fGenerateEnergyOption == "uniform") {
      G4double energy = gRandom->Uniform(fEnergyMin, fEnergyMax);
      if (fEnergyOrMomentum == "energy")
        fParticleGun->SetParticleEnergy(energy);
      else if (fEnergyOrMomentum == "momentum")
        fParticleGun->SetParticleMomentum(energy);
    } else if (fGenerateEnergyOption == "step") {
      G4double energy = (G4double)gRandom->Integer(fEnergyNstep);
      energy = fEnergyMin + energy * fEnergyStepSize;
    }

    // angle
    if (fGenerateAngleOption != "fixed") {
      fPhi = gRandom->Uniform(fPhiMin, fPhiMax);
      // optional fixed phi
      if (fOptionConfiguration.find("fixedPhi") != fOptionConfiguration.end()) {
        fPhi = atof(fOptionConfiguration["fixedPhi"].data()) * deg;
      }
      if (fGenerateAngleOption == "uniformSolid") {
        fDirZ = gRandom->Uniform(cos(fThetaMax), cos(fThetaMin));
        fDirX = sqrt(1.0 - fDirZ * fDirZ) * cos(fPhi);
        fDirY = sqrt(1.0 - fDirZ * fDirZ) * sin(fPhi);
      } else if (fGenerateAngleOption == "uniformTheta") {
        fTheta = gRandom->Uniform(fThetaMin, fThetaMax);
        fDirX = sin(fTheta) * cos(fPhi);
        fDirY = sin(fTheta) * sin(fPhi);
        fDirZ = cos(fTheta);
      } else if (fGenerateAngleOption == "stepTheta") {
        G4double iTheta = (G4double)gRandom->Integer(fThetaNstep);
        fTheta = fThetaMin + iTheta * fThetaStepSize;
        fDirX = sin(fTheta) * cos(fPhi);
        fDirY = sin(fTheta) * sin(fPhi);
        fDirZ = cos(fTheta);
      }
    }
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(fDirX, fDirY, fDirZ));

    gPrimaryParticlePosition = fParticleGun->GetParticlePosition();
    gPrimaryParticleEnergy = fParticleGun->GetParticleEnergy();
    gPrimaryParticleMomentumDirection = fParticleGun->GetParticleMomentumDirection();
    gPrimaryParticlePDG = fParticleGun->GetParticleDefinition()->GetPDGEncoding();
    gPrimaryParticleMass = fParticleGun->GetParticleDefinition()->GetPDGMass();
    fParticleGun->GeneratePrimaryVertex(event);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
