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

#include "KotoEMCalEventAction.hh"
#include "KotoEMCalRunAction.hh"
#include "KotoEMCalEmCalorimeterHit.hh"
#include "KotoEMCalLeadHit.hh"
#include "KotoEMCalCsIHit.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include "KotoEMCalPrimaryGeneratorAction.hh"
#include "TClonesArray.h"
#include "TTree.h"
#include "TObject.h"
#include <Riostream.h>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern G4ThreeVector gPrimaryParticlePosition;
extern G4ThreeVector gPrimaryParticleMomentumDirection;
extern int gPrimaryParticlePDG;
extern double gPrimaryParticleEnergy;
extern double gPrimaryParticleMass;
KotoEMCalEventAction::KotoEMCalEventAction(KotoEMCalRunAction *runAction, TTree *tr)
  : G4UserEventAction(), fRunAction(runAction), fTree(tr)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KotoEMCalEventAction::~KotoEMCalEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KotoEMCalEventAction::BeginOfEventAction(const G4Event*)
{
  SetRunID(fRunAction -> GetRunID());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KotoEMCalEventAction::EndOfEventAction(const G4Event* event)
{
  
  if (fSaveStepLevel){ 
    fEMStepEdep.clear();
    fEMPreStepx.clear();
    fEMPreStepy.clear();
    fEMPreStepz.clear();
    fEMPreStept.clear();
    fEMPostStepx.clear();
    fEMPostStepy.clear();
    fEMPostStepz.clear();
    fEMPostStept.clear();

    fEMParticlePx.clear();
    fEMParticlePy.clear();
    fEMParticlePz.clear();
    fEMParticleTrackID.clear();
    fEMParticleParentID.clear();
    fEMParticleCharge.clear();
    fEMParticleMass.clear();
    fEMParticlePDGID.clear();
  
    
    fLeadStepEdep.clear();
    fLeadPreStepx.clear();
    fLeadPreStepy.clear();
    fLeadPreStepz.clear();
    fLeadPreStept.clear();
    fLeadPostStepx.clear();
    fLeadPostStepy.clear();
    fLeadPostStepz.clear();
    fLeadPostStept.clear();

    fLeadParticlePx.clear();
    fLeadParticlePy.clear();
    fLeadParticlePz.clear();
    fLeadParticleTrackID.clear();
    fLeadParticleParentID.clear();
    fLeadParticleCharge.clear();
    fLeadParticleMass.clear();
    fLeadParticlePDGID.clear();

  }
  EventInfo.eventID = event -> GetEventID();
  if (EventInfo.eventID % 100 == 0 ) G4cout << "Event No. : " << EventInfo.eventID << G4endl;
  
  auto hce = event -> GetHCofThisEvent();
  int iarrayEMHit = 0;
  int iarrayLeadHit = 0;
  int iarrayCsIHit = 0;
  for (int i = 0; i < hce -> GetCapacity(); i++){
    if (hce -> GetHC(i) -> GetSize() == 0) continue;
    G4String iHCName = hce -> GetHC(i) -> GetName();
    // EM Hit
    if (iHCName == "EMCalHitCollection"){
      int nEMHitsInLayer = hce -> GetHC(i) -> GetSize();
      for (int ih = 0; ih < nEMHitsInLayer; ih++){
        auto hit = (KotoEMCalEmCalorimeterHit*) (hce -> GetHC(i) -> GetHit(ih));
        EMHit.cid[iarrayEMHit] = hit -> GetCellID();
        EMHit.lid[iarrayEMHit] = hit -> GetLayerID();
        EMHit.segid[iarrayEMHit] = hit -> GetSegmentID();
        double xx,yy,zz,tt,ee;
        hit -> GetXYZTE(xx,yy,zz,tt,ee);

        EMHit.one[iarrayEMHit] = 1;
        EMHit.x[iarrayEMHit] = xx;
        EMHit.y[iarrayEMHit] = yy;
        EMHit.z[iarrayEMHit] = zz;
        EMHit.t[iarrayEMHit] = tt;
        EMHit.e[iarrayEMHit] = ee;
        if (fSaveStepLevel){
          vector<double> prex, prey, prez, pret;
          vector<double> postx, posty, postz, postt;
          vector<double> stepE;
          vector<double> ppx, ppy, ppz;
          vector<int> trackid, parentid;
          vector<double> charge, mass;
          vector<int> pid;
          hit -> GetPreStepPos(prex,prey,prez,pret);
          hit -> GetPostStepPos(postx,posty,postz,postt);
          hit -> GetStepEdep(stepE);
          hit -> GetParticleTrackInfo(ppx,ppy,ppz,trackid,parentid,charge,mass,pid);
	

          fEMStepEdep.push_back(stepE);
          fEMPreStepx.push_back(prex);
          fEMPreStepy.push_back(prey);
          fEMPreStepz.push_back(prez);
          fEMPreStept.push_back(pret);
          fEMPostStepx.push_back(postx);
          fEMPostStepy.push_back(posty);
          fEMPostStepz.push_back(postz);
          fEMPostStept.push_back(postt);
          /*
            fEMParticlePx.push_back(ppx);
            fEMParticlePy.push_back(ppy);
            fEMParticlePz.push_back(ppz);
            fEMParticleTrackID.push_back(trackid);
            fEMParticleParentID.push_back(parentid);
            fEMParticleCharge.push_back(charge);
            fEMParticleMass.push_back(mass);
            fEMParticlePDGID.push_back(pid);
          */
        }
        iarrayEMHit++;
      }
    }

    // Lead Hit
    else if (iHCName == "LeadHitCollection"){
      auto hit = (KotoEMCalLeadHit*) (hce -> GetHC(i) -> GetHit(0));
      LeadHit.cid[iarrayLeadHit] = hit -> GetCellID();
      LeadHit.lid[iarrayLeadHit] = hit -> GetLayerID();
      double xx,yy,zz,tt,ee;
      hit -> GetXYZTE(xx,yy,zz,tt,ee);
      LeadHit.one[iarrayLeadHit] = 1;
      LeadHit.x[iarrayLeadHit] = xx;
      LeadHit.y[iarrayLeadHit] = yy;
      LeadHit.z[iarrayLeadHit] = zz;
      LeadHit.t[iarrayLeadHit] = tt;
      LeadHit.e[iarrayLeadHit] = ee;	
      
      if (fSaveStepLevel){
        vector<double> prex, prey, prez, pret;
        vector<double> postx, posty, postz, postt;
        vector<double> stepE;
        vector<double> ppx, ppy, ppz;
        vector<int> trackid, parentid;
        vector<double> charge, mass;
        vector<int> pid;
	
        hit -> GetPreStepPos(prex,prey,prez,pret);
        hit -> GetPostStepPos(postx,posty,postz,postt);
        hit -> GetStepEdep(stepE);
        hit -> GetParticleTrackInfo(ppx,ppy,ppz,trackid,parentid,charge,mass,pid);

        fLeadStepEdep.push_back(stepE);
        fLeadPreStepx.push_back(prex);
        fLeadPreStepy.push_back(prey);
        fLeadPreStepz.push_back(prez);
        fLeadPreStept.push_back(pret);
        fLeadPostStepx.push_back(postx);
        fLeadPostStepy.push_back(posty);
        fLeadPostStepz.push_back(postz);
        fLeadPostStept.push_back(postt);
        /*
          fLeadParticlePx.push_back(ppx);
          fLeadParticlePy.push_back(ppy);
          fLeadParticlePz.push_back(ppz);
          fLeadParticleTrackID.push_back(trackid);
          fLeadParticleParentID.push_back(parentid);
          fLeadParticleCharge.push_back(charge);
          fLeadParticleMass.push_back(mass);
          fLeadParticlePDGID.push_back(pid);
        */
      }
      iarrayLeadHit++;
    }

    // CsI Hit
    
    else if (iHCName == "CsIHitCollection"){
      
      auto hit = (KotoEMCalCsIHit*) (hce -> GetHC(i) -> GetHit(0));
      CsIHit.e[iarrayCsIHit] = hit -> GetEdep();
      CsIHit.xid[iarrayCsIHit] = hit -> GetXID();
      CsIHit.yid[iarrayCsIHit] = hit -> GetYID();
      iarrayCsIHit++;
    }
  }
  
  EMHit.nhit = iarrayEMHit;
  LeadHit.nhit = iarrayLeadHit;
  CsIHit.nhit = iarrayCsIHit;
  
  PrimaryParticle.x = gPrimaryParticlePosition.getX();
  PrimaryParticle.y = gPrimaryParticlePosition.getY();
  PrimaryParticle.z = gPrimaryParticlePosition.getZ();
  
  PrimaryParticle.e = gPrimaryParticleEnergy;
  PrimaryParticle.m = gPrimaryParticleMass;
  PrimaryParticle.p = sqrt(PrimaryParticle.e*PrimaryParticle.e - PrimaryParticle.m*PrimaryParticle.m); // magnitude of momentum
  
  PrimaryParticle.px = gPrimaryParticleMomentumDirection.getX() * PrimaryParticle.p;
  PrimaryParticle.py = gPrimaryParticleMomentumDirection.getY() * PrimaryParticle.p;
  PrimaryParticle.pz = gPrimaryParticleMomentumDirection.getZ() * PrimaryParticle.p;
  fTree -> Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KotoEMCalEventAction::SetBranch(){
  if (fSaveStepLevel){
    gInterpreter -> GenerateDictionary("vector<vector<int> >","vector");
    gInterpreter -> GenerateDictionary("vector<vector<double> >","vector");
    //gSystem -> Exec("rm -f AutoDict_vector_vector_*___*");
  }
  fTree -> Branch("eventID",&EventInfo.eventID,"eventID/I");
  fTree -> Branch("runID",&EventInfo.runID,"runID/I");
  fTree -> Branch("randomSeed",&EventInfo.randomSeed,"randomSeed/L");
  
  fTree -> Branch("PrimaryParticle.x",&PrimaryParticle.x,"PrimaryParticle.x/D");
  fTree -> Branch("PrimaryParticle.y",&PrimaryParticle.y,"PrimaryParticle.y/D");
  fTree -> Branch("PrimaryParticle.z",&PrimaryParticle.z,"PrimaryParticle.z/D");

  fTree -> Branch("PrimaryParticle.px",&PrimaryParticle.px,"PrimaryParticle.px/D");
  fTree -> Branch("PrimaryParticle.py",&PrimaryParticle.py,"PrimaryParticle.py/D");
  fTree -> Branch("PrimaryParticle.pz",&PrimaryParticle.pz,"PrimaryParticle.pz/D");
  fTree -> Branch("PrimaryParticle.p",&PrimaryParticle.p,"PrimaryParticle.p/D");
  fTree -> Branch("PrimaryParticle.m",&PrimaryParticle.m,"PrimaryParticle.m/D");
  fTree -> Branch("PrimaryParticle.e",&PrimaryParticle.e,"PrimaryParticle.e/D");
  fTree -> Branch("PrimaryParticle.PDG",&PrimaryParticle.PDG,"PrimaryParticle.PDG/I");

  fTree -> Branch("nEMHit",&EMHit.nhit,"nEMHit/I");
  fTree -> Branch("EMHit.one",EMHit.one,"EMHit.one[nEMHit]/I");
  fTree -> Branch("EMHit.CellID",EMHit.cid,"EMHit.CellID[nEMHit]/I");
  fTree -> Branch("EMHit.LayerID",EMHit.lid,"EMHit.LayerID[nEMHit]/I");
  fTree -> Branch("EMHit.SegmentID",EMHit.segid,"EMHit.SegmentID[nEMHit]/I");
  fTree -> Branch("EMHit.x",EMHit.x,"EMHit.x[nEMHit]/D");
  fTree -> Branch("EMHit.y",EMHit.y,"EMHit.y[nEMHit]/D");
  fTree -> Branch("EMHit.z",EMHit.z,"EMHit.z[nEMHit]/D");
  fTree -> Branch("EMHit.t",EMHit.t,"EMHit.t[nEMHit]/D");
  fTree -> Branch("EMHit.e",EMHit.e,"EMHit.e[nEMHit]/D");
  
  fTree -> Branch("nLeadHit",&LeadHit.nhit,"nLeadHit/I");
  fTree -> Branch("LeadHit.one",LeadHit.one,"LeadHit.one[nLeadHit]/I");
  fTree -> Branch("LeadHit.CellID",LeadHit.cid,"LeadHit.CellID[nLeadHit]/I");
  fTree -> Branch("LeadHit.LayerID",LeadHit.lid,"LeadHit.LayerID[nLeadHit]/I");
  fTree -> Branch("LeadHit.x",LeadHit.x,"LeadHit.x[nLeadHit]/D");
  fTree -> Branch("LeadHit.y",LeadHit.y,"LeadHit.y[nLeadHit]/D");
  fTree -> Branch("LeadHit.z",LeadHit.z,"LeadHit.z[nLeadHit]/D");
  fTree -> Branch("LeadHit.t",LeadHit.t,"LeadHit.t[nLeadHit]/D");
  fTree -> Branch("LeadHit.e",LeadHit.e,"LeadHit.e[nLeadHit]/D");

  fTree -> Branch("nCsIHit",&CsIHit.nhit,"nCsIHit/I");
  fTree -> Branch("CsIHit.xid",CsIHit.xid,"CsIHit.xid[nCsIHit]/I");
  fTree -> Branch("CsIHit.yid",CsIHit.yid,"CsIHit.yid[nCsIHit]/I");
  fTree -> Branch("CsIHit.e",CsIHit.e,"CsIHit.e[nCsIHit]/D");
  
  
  if (fSaveStepLevel){
    fTree -> Branch("EMStepEdep",&fEMStepEdep);
    fTree -> Branch("EMPreStepx",&fEMPreStepx);
    fTree -> Branch("EMPreStepy",&fEMPreStepy);
    fTree -> Branch("EMPreStepz",&fEMPreStepz);
    fTree -> Branch("EMPreStept",&fEMPreStept);
    fTree -> Branch("EMPostStepx",&fEMPostStepx);
    fTree -> Branch("EMPostStepy",&fEMPostStepy);
    fTree -> Branch("EMPostStepz",&fEMPostStepz);
    fTree -> Branch("EMPostStept",&fEMPostStept);
    /*
    fTree -> Branch("EMParticlePx",&fEMParticlePx);
    fTree -> Branch("EMParticlePy",&fEMParticlePy);
    fTree -> Branch("EMParticlePz",&fEMParticlePz);
    fTree -> Branch("EMParticleTrackID",&fEMParticleTrackID);
    fTree -> Branch("EMParticleParentID",&fEMParticleParentID);
    fTree -> Branch("EMParticleCharge",&fEMParticleCharge);
    fTree -> Branch("EMParticleMass",&fEMParticleMass);
    fTree -> Branch("EMParticlePID",&fEMParticlePDGID);
    */
    fTree -> Branch("LeadStepEdep",&fLeadStepEdep);
    fTree -> Branch("LeadPreStepx",&fLeadPreStepx);
    fTree -> Branch("LeadPreStepy",&fLeadPreStepy);
    fTree -> Branch("LeadPreStepz",&fLeadPreStepz);
    fTree -> Branch("LeadPreStept",&fLeadPreStept);
    fTree -> Branch("LeadPostStepx",&fLeadPostStepx);
    fTree -> Branch("LeadPostStepy",&fLeadPostStepy);
    fTree -> Branch("LeadPostStepz",&fLeadPostStepz);
    fTree -> Branch("LeadPostStept",&fLeadPostStept);
    /*
    fTree -> Branch("LeadParticlePx",&fLeadParticlePx);
    fTree -> Branch("LeadParticlePy",&fLeadParticlePy);
    fTree -> Branch("LeadParticlePz",&fLeadParticlePz);
    fTree -> Branch("LeadParticleTrackID",&fLeadParticleTrackID);
    fTree -> Branch("LeadParticleParentID",&fLeadParticleParentID);
    fTree -> Branch("LeadParticleCharge",&fLeadParticleCharge);
    fTree -> Branch("LeadParticleMass",&fLeadParticleMass);
    fTree -> Branch("LeadParticlePID",&fLeadParticlePDGID);
    */ 
  }

}
void KotoEMCalEventAction::SetRunID(G4int RunID){
  EventInfo.runID = RunID;
}
void KotoEMCalEventAction::SetRandomSeed(long seed){
  EventInfo.randomSeed = seed;
}
