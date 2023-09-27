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
/// \file KotoEMCalDetectorConstruction.hh
/// \brief Definition of the KotoEMCalDetectorConstruction class

#ifndef KotoEMCalDetectorConstruction_h
#define KotoEMCalDetectorConstruction_h 1

// c++ std
#include <vector>

// Geant4 class
#include "G4FieldManager.hh"
#include "G4Material.hh"
#include "G4RotationMatrix.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class KotoEMCalMagneticField;

class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4VisAttributes;
class G4GenericMessenger;

/// Detector construction

class KotoEMCalDetectorConstruction : public G4VUserDetectorConstruction {
 public:
  KotoEMCalDetectorConstruction();
  virtual ~KotoEMCalDetectorConstruction();

  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();

 private:
  std::vector<G4VisAttributes*> fVisAttributes;
  // Single Module configuration
  G4int fNLayersInModule;
  G4int fNScintSegmentXYInModule;

  // Module configuration
  G4int fNumberOfModuleX;
  G4int fNumberOfModuleLayer;

  // Absorber
  G4double fAbsorber_size_x;
  G4double fAbsorber_size_y;
  G4double fAbsorber_size_z;

  // Scintillating fiber
  G4double fScintillator_size_x;
  G4double fScintillator_size_y;
  G4double fScintillator_size_z;

  // CsI Box
  G4double fCsI_size_x;
  G4double fCsI_size_y;
  G4double fCsI_size_z;

  G4double fModule_size_x;
  G4double fModule_size_y;
  G4double fModule_size_z;

  G4int fNCsIx;
  G4int fNCsIy;

  G4Material* fMaterial_Absorber;
  G4Material* fMaterial_Scintillator;
  G4Material* fMaterial_CsI;

  G4LogicalVolume* fLogicScint;
  G4LogicalVolume* fLogicAbsorber;
  G4LogicalVolume* fLogicModule;
  G4LogicalVolume* fLogicCsI;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
