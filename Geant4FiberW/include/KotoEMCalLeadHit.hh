#ifndef KotoEMCalLeadHit_h
#define KotoEMCalLeadHit_h 1

// Geant4 class
#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"

class G4AttDef;
class G4AttValue;

/// EM Calorimeter hit
///
/// It records:
/// - the cell ID
/// - the energy deposit
/// - the cell logical volume, its position and rotation

class KotoEMCalLeadHit : public G4VHit {
 public:
  KotoEMCalLeadHit();
  KotoEMCalLeadHit(G4int cellID);
  KotoEMCalLeadHit(const KotoEMCalLeadHit &right);
  virtual ~KotoEMCalLeadHit();

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  G4int GetDetType() const { return fDetType; }  // Lead : 0, Scintillator : 1

  void SetCellID(G4int z) { fCellID = z; }
  G4int GetCellID() const { return fCellID; }

  void SetLayerID(G4int z) { fLayerID = z; }
  G4int GetLayerID() const { return fLayerID; }

  void SetEdep(G4double de) { fEdep = de; }
  void AddEdep(G4double de) { fEdep += de; }
  G4double GetEdep() const { return fEdep; }

  void SetPos(G4ThreeVector xyz) { fPos = xyz; }
  G4ThreeVector GetPos() const { return fPos; }

  void SetStepEdep(std::vector<G4double> stepE) {
    fStepEdep = stepE;
  }

  void GetStepEdep(std::vector<G4double> &stepE) {
    stepE = fStepEdep;
  }

  void SetPreStepPos(std::vector<G4double> prex, std::vector<G4double> prey, std::vector<G4double> prez, std::vector<G4double> pret) {
    fPreStepx = prex;
    fPreStepy = prey;
    fPreStepz = prez;
    fPreStept = pret;
  }

  void SetPostStepPos(std::vector<G4double> postx, std::vector<G4double> posty, std::vector<G4double> postz, std::vector<G4double> postt) {
    fPostStepx = postx;
    fPostStepy = posty;
    fPostStepz = postz;
    fPostStept = postt;
  }

  void SetParticleTrackInfo(std::vector<G4double> ppx, std::vector<G4double> ppy, std::vector<G4double> ppz,
                            std::vector<G4int> trackid, std::vector<G4int> parentid,
                            std::vector<G4double> charge, std::vector<G4double> mass, std::vector<G4int> pid) {
    fParticlePx = ppx;
    fParticlePy = ppy;
    fParticlePz = ppz;
    fParticleTrackID = trackid;
    fParticleParentID = parentid;
    fParticleCharge = charge;
    fParticleMass = mass;
    fParticlePDGID = pid;
  }

  void GetPreStepPos(std::vector<G4double> &prex, std::vector<G4double> &prey, std::vector<G4double> &prez, std::vector<G4double> &pret) {
    prex = fPreStepx;
    prey = fPreStepy;
    prez = fPreStepz;
    pret = fPreStept;
  }

  void GetPostStepPos(std::vector<G4double> &postx, std::vector<G4double> &posty, std::vector<G4double> &postz, std::vector<G4double> &postt) {
    postx = fPostStepx;
    posty = fPostStepy;
    postz = fPostStepz;
    postt = fPostStept;
  }
  void GetParticleTrackInfo(std::vector<G4double> &ppx, std::vector<G4double> &ppy, std::vector<G4double> &ppz,
                            std::vector<G4int> &trackid, std::vector<G4int> &parentid,
                            std::vector<G4double> &charge, std::vector<G4double> &mass, std::vector<G4int> &pid) {
    ppx = fParticlePx;
    ppy = fParticlePy;
    ppz = fParticlePz;
    trackid = fParticleTrackID;
    parentid = fParticleParentID;
    charge = fParticleCharge;
    mass = fParticleMass;
    pid = fParticlePDGID;
  }

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
    G4cout << "(" << fPos.x() << ", " << fPos.y() << ", " << fPos.z() << ", " << fTime << ", " << fEdep << ")" << G4endl;
  }

 private:
  G4int fCellID;
  G4int fLayerID;

  G4double fEdep;
  G4double fTime;
  G4ThreeVector fPos;
  std::vector<G4double> fStepEdep;

  std::vector<G4double> fPreStepx;
  std::vector<G4double> fPreStepy;
  std::vector<G4double> fPreStepz;
  std::vector<G4double> fPreStept;

  std::vector<G4double> fPostStepx;
  std::vector<G4double> fPostStepy;
  std::vector<G4double> fPostStepz;
  std::vector<G4double> fPostStept;

  std::vector<G4double> fParticlePx;
  std::vector<G4double> fParticlePy;
  std::vector<G4double> fParticlePz;
  std::vector<G4int> fParticleTrackID;
  std::vector<G4int> fParticleParentID;
  std::vector<G4double> fParticleCharge;
  std::vector<G4double> fParticleMass;
  std::vector<G4int> fParticlePDGID;

  const G4LogicalVolume *fPLogV;
  const G4int fDetType = 0;  // Lead : 0, Scintillator : 1
};

using KotoEMCalLeadHitsCollection = G4THitsCollection<KotoEMCalLeadHit>;

extern G4ThreadLocal G4Allocator<KotoEMCalLeadHit> *KotoEMCalLeadHitAllocator;

inline void *KotoEMCalLeadHit::operator new(size_t) {
  if (!KotoEMCalLeadHitAllocator) {
    KotoEMCalLeadHitAllocator = new G4Allocator<KotoEMCalLeadHit>;
  }
  return (void *)KotoEMCalLeadHitAllocator->MallocSingle();
}

inline void KotoEMCalLeadHit::operator delete(void *aHit) {
  KotoEMCalLeadHitAllocator->FreeSingle((KotoEMCalLeadHit *)aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
