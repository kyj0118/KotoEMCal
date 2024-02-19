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
#include "KotoEMCalAbsorberSD.hh"
#include "KotoEMCalCsISD.hh"
#include "KotoEMCalDetectorConstruction.hh"
#include "KotoEMCalScintillatorSD.hh"
#include "KotoEMCalTriggerCounterSD.hh"

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

using namespace std;
extern map<string, string> options;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
KotoEMCalDetectorConstruction::KotoEMCalDetectorConstruction()
    : G4VUserDetectorConstruction() {
  // Number of Modules
  fNLayersInModule = 5;
  fNScintSegmentXYInModule = 14;

  fNumberOfModuleLayer = 24;
  fNumberOfModuleX = 16;

  fModuleIntervalXY = 14. * mm;
  fModuleIntervalZ = 6. * mm;

  fNLayersX0Finder = 2;
  fNFibersX0Finder = 16;

  fNTriggerCounter = 2;
  fTrigger_size_x = 10. * cm;
  fTrigger_size_y = 3. * cm;
  fTrigger_size_z = 0.5 * cm;
  fTrigger_interval_z = 50. * cm;

  fNCsIx = 5;
  fNCsIy = 5;
  G4int NumberOfScintillators = fNScintSegmentXYInModule * fNumberOfModuleX;  // Number of scintillator segments in a layer

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
  G4double csi_offset_z = ((G4double)fNumberOfModuleLayer) * fModuleIntervalZ + 0.5 * fCsI_size_z;
  csi_offset_z += 3.0 * cm;
  vector<G4double> vRot;
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

  G4Box* solidTrigger = new G4Box("Trigger",
                                  0.5 * fTrigger_size_x,
                                  0.5 * fTrigger_size_y,
                                  0.5 * fTrigger_size_z);

  G4Box* solidModule = new G4Box("Module",
                                 0.5 * fModule_size_x,
                                 0.5 * fModule_size_y,
                                 0.5 * fModule_size_z);

  G4Box* solidCsI = new G4Box("CsI",
                              0.5 * fCsI_size_x,
                              0.5 * fCsI_size_y,
                              0.5 * fCsI_size_z);

  G4double full_detector_size_z = (csi_offset_z + 0.5 * fCsI_size_z);
  G4Box* solidDetectorMother = new G4Box("DetectorMother",
                                         0.5 * fCsI_size_x * fNCsIx,
                                         0.5 * fCsI_size_y * fNCsIy,
                                         0.5 * full_detector_size_z);

  fLogicScint = new G4LogicalVolume(solidScint, fMaterial_Scintillator, "logicScint");
  fLogicModule = new G4LogicalVolume(solidModule, fMaterial_Scintillator, "logicModule");
  fLogicCsI = new G4LogicalVolume(solidCsI, fMaterial_Scintillator, "logicCsI");
  fLogicAbsorber = new G4LogicalVolume(solidAborber, fMaterial_Absorber, "logicAbsorber");
  fLogicX0Finder = new G4LogicalVolume(solidScint, fMaterial_Scintillator, "LogicX0Finder");
  fLogicTriggerCounter = new G4LogicalVolume(solidTrigger, fMaterial_Scintillator, "LogicTriggerCounter");
  auto fLogicDetectorMother = new G4LogicalVolume(solidDetectorMother, world_mat, "LogicDetectorMother");

  // visualization attributes
  auto visAttributes_Scint = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
  auto visAttributes_Module = new G4VisAttributes(G4Colour(1, 45. / 255., 208. / 255.));  // magenta
  auto visAttributes_Absorber = new G4VisAttributes(G4Colour(0, 0, 0));
  auto visAttributes_CsI = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
  auto visAttributes_MotherDet = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 9));

  visAttributes_MotherDet->SetVisibility(false);
  visAttributes_Module->SetDaughtersInvisible();
  fLogicScint->SetVisAttributes(visAttributes_Scint);
  fLogicAbsorber->SetVisAttributes(visAttributes_Absorber);
  fLogicModule->SetVisAttributes(visAttributes_Module);
  fLogicCsI->SetVisAttributes(visAttributes_CsI);
  fLogicDetectorMother->SetVisAttributes(visAttributes_MotherDet);

  fVisAttributes.push_back(visAttributes_Scint);
  fVisAttributes.push_back(visAttributes_Absorber);
  fVisAttributes.push_back(visAttributes_Module);
  fVisAttributes.push_back(visAttributes_CsI);
  fVisAttributes.push_back(visAttributes_MotherDet);

  // Fill scintillators and absorbers inside the module
  for (int iLayer = 0; iLayer < fNLayersInModule; iLayer++) {
    G4double Sheet_size_z = fScintillator_size_z + fAbsorber_size_z;
    G4double zpos_Absorber = -(0.5 * fModule_size_z) + ((G4double)iLayer) * Sheet_size_z + 0.5 * fAbsorber_size_z;
    G4double zpos_Scint = zpos_Absorber + 0.5 * Sheet_size_z;

    new G4PVPlacement(0,
                      G4ThreeVector(0, 0, zpos_Absorber),
                      fLogicAbsorber,
                      "Absorber",
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

  double rotation_theta = atof(options["detRotationAngle"].c_str()) * deg;
  rotation_theta = -rotation_theta;
  auto RotateDetector = new G4RotationMatrix();
  RotateDetector->rotateY(-rotation_theta);
  // G4ThreeVector RotationArmPosition(0, 0, full_detector_size_z - fCsI_size_z); // rotate axis: CsI surface
  G4ThreeVector RotationArmPosition(0, 0, full_detector_size_z - (fCsI_size_z + 3.0 * cm));  // rear end of EMCal
  G4ThreeVector MotherLogicalPosition(0, 0, 0.5 * full_detector_size_z);
  MotherLogicalPosition = MotherLogicalPosition - RotationArmPosition;
  MotherLogicalPosition.rotateY(rotation_theta);
  MotherLogicalPosition += RotationArmPosition;
  // Place position finding fibers

  G4double pos_z_offset = full_detector_size_z - fCsI_size_z - 33. * cm;
  for (int ii = 0; ii < 2; ii++) {
    G4double X0FinderInterval_z = 50 * cm;
    for (int iLayer = 0; iLayer < fNLayersX0Finder; iLayer++) {
      G4double pos_z = fScintillator_size_z / 2.0 + fScintillator_size_z * ((G4double)(iLayer - fNLayersX0Finder));
      pos_z += pos_z_offset + X0FinderInterval_z * ((G4double)(-ii));
      for (int iScint = 0; iScint < fNFibersX0Finder; iScint++) {
        G4double pos_xy = (((G4double)iScint) - ((G4double)(fNFibersX0Finder - 1)) * 0.5) * fScintillator_size_y;
        G4ThreeVector pos_vector;
        auto rotMatrix = new G4RotationMatrix;
        if (iLayer % 2 == 0) {
          pos_vector = G4ThreeVector(0, pos_xy, pos_z);
          rotMatrix = 0;
        } else {
          pos_vector = G4ThreeVector(pos_xy, 0, pos_z);
          rotMatrix = pRot;
        }
        // rotMatrix->rotateY(-rotation_theta);
        new G4PVPlacement(rotMatrix,
                          pos_vector,
                          fLogicX0Finder,
                          "X0Finder",
                          logicWorld,
                          false,
                          iLayer * fNFibersX0Finder + iScint + fNFibersX0Finder * fNLayersX0Finder * ii,
                          false);
      }
    }
  }
  // Place trigger counters
  for (int itrig = 0; itrig < fNTriggerCounter; itrig++) {
    G4double triggerCounterOffset_z = -50 * cm;
    G4double pos_z = ((G4double)(itrig - fNTriggerCounter)) * fTrigger_interval_z + triggerCounterOffset_z;
    G4ThreeVector OriginalPosition(0, 0, pos_z);
    // G4ThreeVector RevisedPosition = OriginalPosition.rotateY(rotation_theta);
    new G4PVPlacement(0,
                      OriginalPosition,
                      fLogicTriggerCounter,
                      "Trigger",
                      logicWorld,
                      false,
                      itrig,
                      false);
  }

  // Place detector mother for rotation
  new G4PVPlacement(RotateDetector,
                    MotherLogicalPosition,
                    fLogicDetectorMother,
                    "DetectorMother",
                    logicWorld,
                    false,
                    0,
                    false);

  // Place modules inside the detector mother
  for (int iLayer = 0; iLayer < fNumberOfModuleLayer; iLayer++) {
    G4double module_pos_z = fModule_size_z / 2.0 + fModuleIntervalZ * ((G4double)iLayer);
    module_pos_z -= 0.5 * full_detector_size_z;
    for (int imod = 0; imod < fNumberOfModuleX; imod++) {
      G4double module_pos_xy = (((G4double)imod) - ((G4double)(fNumberOfModuleX - 1)) * 0.5) * fModuleIntervalXY;
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
                        fLogicDetectorMother,
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
                        G4ThreeVector(csi_offset_x, csi_offset_y, csi_offset_z - full_detector_size_z * 0.5),
                        fLogicCsI,
                        "CsI",
                        fLogicDetectorMother,
                        false,
                        i * fNCsIy + j,
                        false);
    }
  }

  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void KotoEMCalDetectorConstruction::ConstructSDandField() {
  // Sensitive Detector

  auto ScintillatorSD = new KotoEMCalScintillatorSD("ScintSD", fNumberOfModuleLayer, fNumberOfModuleX);
  G4SDManager::GetSDMpointer()->AddNewDetector(ScintillatorSD);
  fLogicScint->SetSensitiveDetector(ScintillatorSD);

  auto AbsorberSD = new KotoEMCalAbsorberSD("AbsorberSD", fNumberOfModuleLayer, fNumberOfModuleX);
  G4SDManager::GetSDMpointer()->AddNewDetector(AbsorberSD);
  fLogicAbsorber->SetSensitiveDetector(AbsorberSD);

  auto CsISD = new KotoEMCalCsISD("CsISD", fNCsIx, fNCsIy);
  G4SDManager::GetSDMpointer()->AddNewDetector(CsISD);
  fLogicCsI->SetSensitiveDetector(CsISD);

  auto X0FinderSD = new KotoEMCalTriggerCounterSD("X0FinderSD", fNLayersX0Finder, fNFibersX0Finder * 2);
  G4SDManager::GetSDMpointer()->AddNewDetector(X0FinderSD);
  fLogicX0Finder->SetSensitiveDetector(X0FinderSD);

  auto TriggerCounterSD = new KotoEMCalTriggerCounterSD("TriggerCounterSD", fNTriggerCounter, 1);
  G4SDManager::GetSDMpointer()->AddNewDetector(TriggerCounterSD);
  fLogicTriggerCounter->SetSensitiveDetector(TriggerCounterSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
