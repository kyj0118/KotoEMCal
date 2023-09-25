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
/// \file KotoEMCalEmCalorimeterSD.cc
/// \brief Implementation of the KotoEMCalEmCalorimeterSD class

#include "KotoEMCalLeadSD.hh"
#include "KotoEMCalLeadHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern bool gSaveStepLevel;

KotoEMCalLeadSD::KotoEMCalLeadSD(G4String name, G4int layerNumber)
  : G4VSensitiveDetector(name), fNameSD(name), fLayerId(layerNumber), fHitsCollection(nullptr), fHCID(-1)
{
  fEdep = 0.0; fEweightedx = 0.0; fEweightedy = 0.0; fEweightedz = 0.0; fEweightedt = 0.0;
  collectionName.insert("LeadHitCollection"); 

  fStepEdep.clear();

  fPreStepx.clear();
  fPreStepy.clear();
  fPreStepz.clear();
  fPreStept.clear();
  
  fPostStepx.clear();
  fPostStepy.clear();
  fPostStepz.clear();
  fPostStept.clear();

  fParticlePx.clear();
  fParticlePy.clear();
  fParticlePz.clear();
  fParticleTrackID.clear();
  fParticleParentID.clear();
  fParticleCharge.clear();
  fParticleMass.clear();
  fParticlePDGID.clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KotoEMCalLeadSD::~KotoEMCalLeadSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KotoEMCalLeadSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection = new KotoEMCalLeadHitsCollection(fNameSD,collectionName[0]);
  if (fHCID<0) {
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }
  hce->AddHitsCollection(fHCID,fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool KotoEMCalLeadSD::ProcessHits(G4Step*step, G4TouchableHistory*)
{
  auto edep = step->GetTotalEnergyDeposit();
  if (edep != 0.0){
    auto prepoint = step -> GetPreStepPoint();
    auto postpoint = step -> GetPostStepPoint();
    
    G4double prex = (prepoint -> GetPosition()).x();
    G4double prey = (prepoint -> GetPosition()).y();
    G4double prez = (prepoint -> GetPosition()).z();
    G4double pret = prepoint -> GetGlobalTime();
    
    G4double postx = (postpoint -> GetPosition()).x();
    G4double posty = (postpoint -> GetPosition()).y();
    G4double postz = (postpoint -> GetPosition()).z();
    G4double postt = prepoint -> GetGlobalTime();
    
    G4double x = (prex + postx)/2.0;
    G4double y = (prey + posty)/2.0;
    G4double z = (prez + postz)/2.0;
    G4double t = (pret + postt)/2.0;

    fEdep += edep;
    fEweightedx += x * edep;
    fEweightedy += y * edep;
    fEweightedz += z * edep;
    fEweightedt += t * edep; 

    if (gSaveStepLevel == true){    
      fStepEdep.push_back(edep);
    
      fPreStepx.push_back(prex);
      fPreStepy.push_back(prey);
      fPreStepz.push_back(prez);
      fPreStept.push_back(pret);

      fPostStepx.push_back(postx);
      fPostStepy.push_back(posty);
      fPostStepz.push_back(postz);
      fPostStept.push_back(postt);

      // particle info
      G4Track *tr = step -> GetTrack();
      const G4ParticleDefinition *pdef = tr -> GetParticleDefinition();
      G4ThreeVector pp = tr -> GetMomentum();
    
      G4double ppx = pp.x();
      G4double ppy = pp.y();
      G4double ppz = pp.z();
      G4int trackid = tr -> GetTrackID();
      G4int parentid = tr -> GetParentID();
      G4double pcharge = pdef -> GetPDGCharge();
      G4double pmass = pdef -> GetPDGMass();
      G4int pid = pdef -> GetPDGEncoding();

      fParticlePx.push_back(ppx);
      fParticlePy.push_back(ppy);
      fParticlePz.push_back(ppz);
      fParticleTrackID.push_back(trackid);
      fParticleParentID.push_back(parentid);
      fParticleCharge.push_back(pcharge);
      fParticleMass.push_back(pmass);
      fParticlePDGID.push_back(pid);
    }
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KotoEMCalLeadSD::EndOfEvent(G4HCofThisEvent* hce){
  if (fEdep == 0) return;
  else {
    fEweightedx/=fEdep;
    fEweightedy/=fEdep;
    fEweightedz/=fEdep;
    fEweightedt/=fEdep;
    
    fHitsCollection->insert(new KotoEMCalLeadHit(fHCID));
    auto hit = (KotoEMCalLeadHit*) ((hce -> GetHC(fHCID)) -> GetHit(0));
    hit -> SetXYZTE(fEweightedx, fEweightedy, fEweightedz, fEweightedt, fEdep);
    hit -> SetLayerID(fLayerId);

    if (gSaveStepLevel == true){
      hit -> SetPreStepPos(fPreStepx,fPreStepy,fPreStepz,fPreStept);
      hit -> SetPostStepPos(fPostStepx,fPostStepy,fPostStepz,fPostStept);
      hit -> SetStepEdep(fStepEdep);
      
      hit -> SetParticleTrackInfo(fParticlePx,fParticlePy,fParticlePz,fParticleTrackID,fParticleParentID,
				  fParticleCharge,fParticleMass,fParticlePDGID);

      fStepEdep.clear();
    
      fPreStepx.clear();
      fPreStepy.clear();
      fPreStepz.clear();
      fPreStept.clear();
    
      fPostStepx.clear();
      fPostStepy.clear();
      fPostStepz.clear();
      fPostStept.clear();

      fParticlePx.clear();
      fParticlePy.clear();
      fParticlePz.clear();
      fParticleTrackID.clear();
      fParticleParentID.clear();
      fParticleCharge.clear();
      fParticleMass.clear();
      fParticlePDGID.clear();
    }
    fEdep = 0.0; fEweightedx = 0.0; fEweightedy = 0.0; fEweightedz = 0.0; fEweightedt = 0.0;
  }
}
