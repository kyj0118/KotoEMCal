#ifndef KotoEMCalAbsorberHit_h
#define KotoEMCalAbsorberHit_h 1

// Geant4 class
#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"

class G4AttDef;
class G4AttValue;

class KotoEMCalAbsorberHit : public G4VHit {
 public:
  KotoEMCalAbsorberHit();
  KotoEMCalAbsorberHit(G4int moduleID);
  KotoEMCalAbsorberHit(const KotoEMCalAbsorberHit &right);
  virtual ~KotoEMCalAbsorberHit();

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

using KotoEMCalAbsorberHitsCollection = G4THitsCollection<KotoEMCalAbsorberHit>;

extern G4ThreadLocal G4Allocator<KotoEMCalAbsorberHit> *KotoEMCalAbsorberHitAllocator;

inline void *KotoEMCalAbsorberHit::operator new(size_t) {
  if (!KotoEMCalAbsorberHitAllocator) {
    KotoEMCalAbsorberHitAllocator = new G4Allocator<KotoEMCalAbsorberHit>;
  }
  return (void *)KotoEMCalAbsorberHitAllocator->MallocSingle();
}

inline void KotoEMCalAbsorberHit::operator delete(void *aHit) {
  KotoEMCalAbsorberHitAllocator->FreeSingle((KotoEMCalAbsorberHit *)aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif