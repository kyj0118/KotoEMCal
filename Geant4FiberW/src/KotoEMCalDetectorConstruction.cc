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
/// \file KotoEMCalDetectorConstruction.cc
/// \brief Implementation of the KotoEMCalDetectorConstruction class

// Root class
#include "TString.h"

// This project class
#include "KotoEMCalCsISD.hh"
#include "KotoEMCalDetectorConstruction.hh"
#include "KotoEMCalEmCalorimeterSD.hh"
#include "KotoEMCalLeadSD.hh"

// Geant4 class
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Element.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4UserLimits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"
#include "G4tgbRotationMatrix.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int gNumberOfScintillators;
KotoEMCalDetectorConstruction::KotoEMCalDetectorConstruction()
    : G4VUserDetectorConstruction() {
  // Number of Modules
  fNLayersInModule = 5;
  fNScintSegmentXYInModule = 14;

  fNumberOfModuleLayer = 24;
  fNumberOfModuleX = 16;

  fNCsIx = 5;
  fNCsIy = 5;
  G4int NumberOfScintillators = fNScintSegmentXYInModule * fNumberOfModuleX;  // Number of scintillator segments in a layer
  gNumberOfScintillators = NumberOfScintillators;

  // Size
  fAbsorber_size_x = 25. * cm;

  fAbsorber_size_y = ((G4double)fNScintSegmentXYInModule) * mm;
  fAbsorber_size_z = 0.15 * mm;

  fScintillator_size_x = fAbsorber_size_x;
  fScintillator_size_y = 1. * mm;
  fScintillator_size_z = 1. * mm;

  fCsI_size_x = 7.0 * cm;
  fCsI_size_y = 7.0 * cm;
  fCsI_size_z = 30.0 * cm;

  fModule_size_x = fAbsorber_size_x;
  fModule_size_y = fAbsorber_size_y;
  fModule_size_z = (fAbsorber_size_z + fScintillator_size_z) * ((G4double)fNLayersInModule);

  // Materials
  auto nist = G4NistManager::Instance();

  fMaterial_Absorber = nist->FindOrBuildMaterial("G4_W");
  fMaterial_Scintillator = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  fMaterial_CsI = nist->FindOrBuildMaterial("G4_CESIUM_IODIDE");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KotoEMCalDetectorConstruction::~KotoEMCalDetectorConstruction() {
  for (auto visAttributes : fVisAttributes) {
    delete visAttributes;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* KotoEMCalDetectorConstruction::Construct() {
  G4NistManager* nist = G4NistManager::Instance();

  // -----------------------------------------------------
  // World

  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4double world_size = 5. * m;

  G4Box* solidWorld =
      new G4Box("World",            // its name
                0.5 * world_size,   // half x
                0.5 * world_size,   // half y
                0.5 * world_size);  // half z

  G4LogicalVolume* logicWorld =
      new G4LogicalVolume(solidWorld,  // its solid
                          world_mat,   // its material
                          "World");    // its name

  G4VPhysicalVolume* physWorld =
      new G4PVPlacement(0,                // no rotation
                        G4ThreeVector(),  // at (0,0,0)
                        logicWorld,       // its logical volume
                        "World",          // its name
                        0,                // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        false);           // overlaps checking

  // -----------------------------------------------------
  G4double layer_gap_z = 6.0 * mm;
  G4double csi_offset_z = ((G4double)fNumberOfModuleLayer) * layer_gap_z + 0.5 * fCsI_size_z;

  std::vector<G4double> vRot;
  vRot.clear();
  vRot.push_back(0.);
  vRot.push_back(1.);
  vRot.push_back(0.);
  vRot.push_back(-1.);
  vRot.push_back(0.);
  vRot.push_back(0.);
  vRot.push_back(0.);
  vRot.push_back(0.);
  vRot.push_back(1.);

  G4tgbRotationMatrix* RotBuilder = new G4tgbRotationMatrix();
  G4RotationMatrix* pRot = (G4RotationMatrix*)RotBuilder->BuildG4RotMatrixFrom9(vRot);

  G4Box* solidAborber = new G4Box("Absorber",
                                  0.5 * fAbsorber_size_x,
                                  0.5 * fAbsorber_size_y,
                                  0.5 * fAbsorber_size_z);

  G4Box* solidScint = new G4Box("Scint",
                                0.5 * fScintillator_size_x,
                                0.5 * fScintillator_size_y,
                                0.5 * fScintillator_size_z);

  G4Box* solidModule = new G4Box("Module",
                                 0.5 * fModule_size_x,
                                 0.5 * fModule_size_y,
                                 0.5 * fModule_size_z);

  G4Box* solidCsI = new G4Box("CsI",
                              0.5 * fCsI_size_x,
                              0.5 * fCsI_size_y,
                              0.5 * fCsI_size_z);

  fLogicScint = new G4LogicalVolume(solidScint, fMaterial_Scintillator, "logicScint");
  fLogicAbsorber = new G4LogicalVolume(solidAborber, fMaterial_Absorber, "logicAbsorber");
  fLogicModule = new G4LogicalVolume(solidModule, fMaterial_Scintillator, "logicModule");
  fLogicCsI = new G4LogicalVolume(solidCsI, fMaterial_Scintillator, "logicCsI");

  // visualization attributes
  auto visAttributes_Scint = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
  auto visAttributes_Module = new G4VisAttributes(G4Colour(1, 45. / 255., 208. / 255.));  // magenta
  auto visAttributes_Absorber = new G4VisAttributes(G4Colour(0, 0, 0));
  auto visAttributes_CsI = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));

  visAttributes_Module->SetDaughtersInvisible();
  // visAttributes_Module -> SetColor(G4Colour(255, 45, 208));
  fLogicScint->SetVisAttributes(visAttributes_Scint);
  fLogicAbsorber->SetVisAttributes(visAttributes_Absorber);
  fLogicModule->SetVisAttributes(visAttributes_Module);
  fLogicCsI->SetVisAttributes(visAttributes_CsI);

  fVisAttributes.push_back(visAttributes_Scint);
  fVisAttributes.push_back(visAttributes_Absorber);
  fVisAttributes.push_back(visAttributes_Module);
  fVisAttributes.push_back(visAttributes_CsI);

  // Fill scintillators and absorbers inside the module
  for (int iLayer = 0; iLayer < fNLayersInModule; iLayer++) {
    G4double Sheet_size_z = fScintillator_size_z + fAbsorber_size_z;
    G4double zpos_Absorber = -(0.5 * fModule_size_z) + ((G4double)iLayer) * Sheet_size_z + 0.5 * fAbsorber_size_z;
    G4double zpos_Scint = zpos_Absorber + 0.5 * Sheet_size_z;

    new G4PVPlacement(0,
                      G4ThreeVector(0, 0, zpos_Absorber),
                      fLogicAbsorber,
                      "Scint",
                      fLogicModule,
                      false,
                      iLayer,  // copy No.
                      false);

    for (int iScint = 0; iScint < fNScintSegmentXYInModule; iScint++) {
      G4double xpos = 0;
      G4double ypos = (((G4double)iScint) - ((G4double)(fNScintSegmentXYInModule - 1)) * 0.5) * fScintillator_size_y;

      new G4PVPlacement(0,
                        G4ThreeVector(xpos, ypos, zpos_Scint),
                        fLogicScint,
                        "Scint",
                        fLogicModule,
                        false,
                        iScint + fNScintSegmentXYInModule * iLayer,  // copy No.
                        false);
    }
  }

  // Place modules on the world
  for (int iLayer = 0; iLayer < fNumberOfModuleLayer; iLayer++) {
    G4double module_pos_z = fModule_size_z / 2.0 + layer_gap_z * ((G4double)iLayer);
    for (int imod = 0; imod < fNumberOfModuleX; imod++) {
      G4double module_pos_xy = (((G4double)imod) - ((G4double)(fNumberOfModuleX - 1)) * 0.5) * fModule_size_y;
      G4ThreeVector pos_vector;
      auto rotMatrix = new G4RotationMatrix;
      if (iLayer % 2 == 0) {
        pos_vector = G4ThreeVector(0, module_pos_xy, module_pos_z);
        rotMatrix = 0;
      } else {
        pos_vector = G4ThreeVector(module_pos_xy, 0, module_pos_z);
        rotMatrix = pRot;
      }
      new G4PVPlacement(rotMatrix,
                        pos_vector,
                        fLogicModule,
                        "Module",
                        logicWorld,
                        false,
                        iLayer * fNumberOfModuleX + imod,
                        false);
    }
  }

  for (int i = 0; i < fNCsIx; i++) {
    G4double csi_offset_x = (((G4double)i) - ((G4double)(fNCsIx - 1)) / 2.0) * fCsI_size_x;
    for (int j = 0; j < fNCsIy; j++) {
      G4double csi_offset_y = (((G4double)j) - ((G4double)(fNCsIy - 1)) / 2.0) * fCsI_size_y;
      new G4PVPlacement(0,
                        G4ThreeVector(csi_offset_x, csi_offset_y, csi_offset_z),
                        fLogicCsI,
                        "CsI",
                        logicWorld,
                        false,
                        0,
                        false);
    }
  }

  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void KotoEMCalDetectorConstruction::ConstructSDandField() {
  // Sensitive Detector
  auto ScintillatorSD = new KotoEMCalEmCalorimeterSD("ScintSD", 0, 0);
  G4SDManager::GetSDMpointer()->AddNewDetector(ScintillatorSD);
  fLogicScint->SetSensitiveDetector(ScintillatorSD);

  auto AbsorberSD = new KotoEMCalLeadSD("AbsorberSD", 0);
  G4SDManager::GetSDMpointer()->AddNewDetector(AbsorberSD);
  fLogicAbsorber->SetSensitiveDetector(AbsorberSD);

  auto CsISD = new KotoEMCalCsISD("CsISD", 0, 0);
  G4SDManager::GetSDMpointer()->AddNewDetector(AbsorberSD);
  fLogicCsI->SetSensitiveDetector(AbsorberSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
