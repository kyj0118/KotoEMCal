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

#include "KotoEMCalDetectorConstruction.hh"
#include "KotoEMCalEmCalorimeterSD.hh"
#include "KotoEMCalLeadSD.hh"
#include "KotoEMCalCsISD.hh"

#include "G4TransportationManager.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"
#include "G4tgbRotationMatrix.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

// Root classes
#include "TString.h"
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int gNumberOfScintillators;
KotoEMCalDetectorConstruction::KotoEMCalDetectorConstruction()
  : G4VUserDetectorConstruction()
{
  NumberOfLayers = 120;
  NumberOfScintillators = 250; // Number of scintillator segments in a layer
  NumberOfLayersGrouping = 5;
  gNumberOfScintillators = NumberOfScintillators;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KotoEMCalDetectorConstruction::~KotoEMCalDetectorConstruction()
{
  
  for (auto visAttributes: fVisAttributes) {
    delete visAttributes;
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* KotoEMCalDetectorConstruction::Construct(){
  G4NistManager* nist = G4NistManager::Instance();
  
  // -----------------------------------------------------
  // World

  G4Material* world_mat = nist -> FindOrBuildMaterial("G4_AIR");
  G4double world_size = 5.*m;
  
  G4Box* solidWorld =    
    new G4Box("World",                       // its name
              0.5*world_size,                // half x
              0.5*world_size,                // half y
              0.5*world_size);               // half z
  
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      false);                //overlaps checking


  // -----------------------------------------------------

  // Detector

  // Materials
  G4Material* Material_W = nist -> FindOrBuildMaterial("G4_W");
  G4Material* Material_Scint = nist -> FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  G4Material* Material_CsI = nist -> FindOrBuildMaterial("G4_CESIUM_IODIDE");
  
  // Tungsten plate
  //G4double W_size_x = 28.5*cm;
  G4double W_size_x = 25*cm;
  //G4double W_size_x = 1.*cm;
  G4double W_size_y = W_size_x;
  G4double W_size_z = 0.15 * mm; 
  //G4double W_size_z = 0.17 * mm; 
  
  // EJ 200 Scintillator
  G4double scint_size_x = W_size_x;
  G4double scint_size_y = 1.*mm;
  G4double scint_size_z = 1.*mm;

  // CsI Box
  G4double csi_size_x = 7.0*cm;
  G4double csi_size_y = 7.0*cm;
  G4double csi_size_z = 30.0 *cm; // ~16X0 length

  // Number of csi cells nx x ny array
  G4int nx_csi = 5;
  G4int ny_csi = 5;
  
  G4double layer_gap_z = W_size_z + scint_size_z;
  
  G4double csi_offset_z = ((G4double) NumberOfLayers)*layer_gap_z + 0.5 * csi_size_z;
  
  
  
  std::vector<G4double> vRot;
  vRot.clear();
  vRot.push_back(0.);vRot.push_back(1.);vRot.push_back(0.);
  vRot.push_back(-1.);vRot.push_back(0.);vRot.push_back(0.);
  vRot.push_back(0.);vRot.push_back(0.);vRot.push_back(1.);
  
  G4tgbRotationMatrix* RotBuilder = new G4tgbRotationMatrix();
  G4RotationMatrix* pRot = (G4RotationMatrix*) RotBuilder->BuildG4RotMatrixFrom9(vRot); 
  
  
  G4Box* solidW = new G4Box("W",
			    0.5*W_size_x,
			    0.5*W_size_y,
			    0.5*W_size_z);
  
  G4Box* solidScint = new G4Box("Scint",
				0.5*scint_size_x,
				0.5*scint_size_y,
				0.5*scint_size_z);

  G4Box* solidScintPlate = new G4Box("ScintPlate",
				     0.5*scint_size_x,
				     0.5*scint_size_x,
				     0.5*scint_size_z);
  
  G4Box* solidCsI = new G4Box("CsI",
			    0.5*csi_size_x,
			    0.5*csi_size_y,
			    0.5*csi_size_z);

  
  auto visAttributes_Scint = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  
  
  for (int iLayer = 0; iLayer < NumberOfLayers; iLayer++){
    G4double dLayer = (G4double) iLayer;
    
    G4double W_offset_z = 0.5 * W_size_z + layer_gap_z * dLayer;
    G4double scint_offset_z = W_offset_z + 0.5 * layer_gap_z;
    
    TString str_leadname_tmp = Form("W%d",iLayer);
    G4String str_leadname = str_leadname_tmp.Data();
    
    TString str_scintname_tmp = Form("Scint%d",iLayer);
    G4String str_scintname = str_scintname_tmp.Data();
    
    G4LogicalVolume* logicScintPlate =  new G4LogicalVolume(solidScintPlate, Material_Scint, str_scintname);
    G4LogicalVolume* logicScint =  new G4LogicalVolume(solidScint, Material_Scint, str_scintname); 

    G4int iLayerGroup = iLayer/NumberOfLayersGrouping;
    if (iLayerGroup % 2 == 0){
      auto pPhysScint = new G4PVPlacement(0,
					  G4ThreeVector(0,0,scint_offset_z),
					  logicScintPlate,
					  str_scintname,
					  logicWorld,
					  false,
					  iLayer,
					  false);
      G4VPhysicalVolume* div_ = new G4PVReplica("div_", logicScint, pPhysScint, kYAxis, NumberOfScintillators, scint_size_y);
    }
    else {
      auto pPhysScint = new G4PVPlacement(pRot,
					  G4ThreeVector(0,0,scint_offset_z),
					  logicScintPlate,
					  str_scintname,
					  logicWorld,
					  false,
					  iLayer,
					  false);
      G4VPhysicalVolume* div_ = new G4PVReplica("div_", logicScint, pPhysScint, kYAxis, NumberOfScintillators, scint_size_y);
    }
    
    
    G4LogicalVolume* logicW = new G4LogicalVolume(solidW, Material_W, str_leadname);
    new G4PVPlacement(0,
		      G4ThreeVector(0,0,W_offset_z),
		      logicW,
		      str_leadname,
		      logicWorld,
		      false,
		      0,
		      false);
    
    
    // visualization attributes 
    logicScintPlate->SetVisAttributes(visAttributes_Scint);
    fVisAttributes.push_back(visAttributes_Scint);  
  }
  
  for (int i = 0; i < nx_csi; i++){
    G4double csi_offset_x = (((G4double)i) - ((G4double)(nx_csi-1))/2.0)*csi_size_x;
    for (int j = 0; j < ny_csi; j++){
      G4double csi_offset_y = (((G4double)j) - ((G4double)(ny_csi-1))/2.0)*csi_size_y;

      G4String str_csiname = Form("CsI_%d_%d",i,j);

      G4LogicalVolume* logicCsI = new G4LogicalVolume(solidCsI, Material_CsI, str_csiname);
      new G4PVPlacement(0,
			G4ThreeVector(csi_offset_x,csi_offset_y,csi_offset_z),
			logicCsI,
			str_csiname,
			logicWorld,
			false,
			0,
			false);
    }
  }
  return physWorld;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void KotoEMCalDetectorConstruction::ConstructSDandField()
{
  //Sensitive Detector
  for (int iLayer = 0 ; iLayer < NumberOfLayers; iLayer++){
    TString str_leadname_tmp = Form("W%d",iLayer);
    G4String str_leadname = str_leadname_tmp.Data();
    G4String str_SDnameLead = (str_leadname_tmp+"_SD").Data();
    
    TString str_scintname_tmp = Form("Scint%d",iLayer);
    G4String str_scintname = str_scintname_tmp.Data();
    G4String str_SDnameScint = (str_scintname_tmp+"_SD").Data();
    
    KotoEMCalEmCalorimeterSD* ScintillatorSD = new KotoEMCalEmCalorimeterSD(str_SDnameScint,iLayer,0);
    G4SDManager::GetSDMpointer() -> AddNewDetector(ScintillatorSD);
    SetSensitiveDetector(str_scintname, ScintillatorSD, true);
    
    // W
    KotoEMCalLeadSD* LeadSD = new KotoEMCalLeadSD(str_SDnameLead,iLayer);
    G4SDManager::GetSDMpointer() -> AddNewDetector(LeadSD);
    SetSensitiveDetector(str_leadname, LeadSD, true);
  }

  // CsI
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 5; j++){
      G4String str_csiname = Form("CsI_%d_%d",i,j);
      G4String str_csinameSD = Form("CsISD_%d_%d",i,j);
      
      KotoEMCalCsISD* CsISD = new KotoEMCalCsISD(str_csinameSD,i,j);
      G4SDManager::GetSDMpointer() -> AddNewDetector(CsISD);
      SetSensitiveDetector(str_csiname, CsISD, true);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

