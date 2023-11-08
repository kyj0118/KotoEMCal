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
/// \file KotoEMCalScintillatorHit.hh
/// \brief Definition of the KotoEMCalScintillatorHit class

#ifndef KotoEMCalScintillatorHit_h
#define KotoEMCalScintillatorHit_h 1

// Geant4 class
#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"

class G4AttDef;
class G4AttValue;

class KotoEMCalScintillatorHit : public G4VHit {
 public:
  KotoEMCalScintillatorHit();
  KotoEMCalScintillatorHit(G4int moduleID);
  KotoEMCalScintillatorHit(const KotoEMCalScintillatorHit &right);
  virtual ~KotoEMCalScintillatorHit();

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  void SetModuleID(G4int z) { fModuleID = z; }
  G4int GetModuleID() const { return fModuleID; }

  void SetLayerID(G4int z) { fLayerID = z; }
  G4int GetLayerID() const { return fLayerID; }

  void SetSegmentID(G4int z) { fSegmentID = z; }
  G4int GetSegmentID() const { return fSegmentID; }

  void SetEdep(G4double de) { fEdep = de; }
  void AddEdep(G4double de) { fEdep += de; }
  G4double GetEdep() const { return fEdep; }

  void SetPos(G4ThreeVector xyz) { fPos = xyz; }
  G4ThreeVector GetPos() const { return fPos; }

  void SetLogV(G4LogicalVolume *val) { fPLogV = val; }
  const G4LogicalVolume *GetLogV() const { return fPLogV; }

  void SetXYZTE(G4double x, G4double y, G4double z, G4double t, G4double e) {
    fPos.setX(x);
    fPos.setY(y);
    fPos.setZ(z);
    fTime = t;
    fEdep = e;
  };

  void GetXYZTE(G4double &x, G4double &y, G4double &z, G4double &t, G4double &e) const {
    x = fPos.x();
    y = fPos.y();
    z = fPos.z();
    t = fTime;
    e = fEdep;
  };
  void Print() {
    std::cout << "(" << fPos.x() << ", " << fPos.y() << ", " << fPos.z() << ", " << fTime << ", " << fEdep << ")" << std::endl;
  }

 private:
  G4int fModuleID;
  G4int fLayerID;
  G4int fSegmentID;

  G4double fEdep;
  G4double fTime;
  G4ThreeVector fPos;

  const G4LogicalVolume *fPLogV;
};

using KotoEMCalScintillatorHitsCollection = G4THitsCollection<KotoEMCalScintillatorHit>;

extern G4ThreadLocal G4Allocator<KotoEMCalScintillatorHit> *KotoEMCalScintillatorHitAllocator;

inline void *KotoEMCalScintillatorHit::operator new(size_t) {
  if (!KotoEMCalScintillatorHitAllocator) {
    KotoEMCalScintillatorHitAllocator = new G4Allocator<KotoEMCalScintillatorHit>;
  }
  return (void *)KotoEMCalScintillatorHitAllocator->MallocSingle();
}

inline void KotoEMCalScintillatorHit::operator delete(void *aHit) {
  KotoEMCalScintillatorHitAllocator->FreeSingle((KotoEMCalScintillatorHit *)aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
