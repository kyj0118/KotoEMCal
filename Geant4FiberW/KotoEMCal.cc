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
/// \file exampleKotoEMCal.cc
/// \brief Main program of the analysis/KotoEMCal example

// Root class
#include "TFile.h"
#include "TObjString.h"
#include "TRandom3.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

// This project class
#include "KotoEMCalActionInitialization.hh"
#include "KotoEMCalDetectorConstruction.hh"
#include "KotoEMCalPrimaryGeneratorAction.hh"

// Genat4 class
#include "FTFP_BERT.hh"
#include "G4RunManager.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"
#include "globals.hh"

// c++ std
#include "time.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;
map<string, string> options;
map<string, string> parseOptions(int argc, char** argv);

int main(int argc, char** argv) {
  // Option parsing
  options = parseOptions(argc, argv);
  if (argc == 1) {
    options["useGPS"] = "true";
    options["mac"] = "vis.mac";
    options["seed"] = "1";
    options["notSaveTree"] = "true";
    options["detRotationAngle"] = "30";
  }

  for (auto [key, value] : options) {
    cout << key << " : " << value << endl;
  }

  TString strOutFile = options["outFileName"].c_str();
  // Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
  long randomSeed = atol(options["seed"].c_str());
  CLHEP::HepRandom::setTheSeed(randomSeed);
  gRandom->SetSeed(randomSeed);

  if (strOutFile == "") strOutFile = "DefaultOutputName.root";
  if (!strOutFile.EndsWith(".root")) {
    strOutFile += ".root";
  }

  // if (options["visualizationOnly"] != "true") {
  auto tf = new TFile(strOutFile, "RECREATE");
  auto tr = new TTree("tree", "Geant4 output");
  tr->SetAutoSave();
  auto physicsList = new FTFP_BERT(0);
  G4RunManager* runManager = new G4RunManager;
  runManager->SetUserInitialization(physicsList);
  runManager->SetUserInitialization(new KotoEMCalActionInitialization(tr));
  runManager->SetUserInitialization(new KotoEMCalDetectorConstruction());
  auto primaryGenerator = new KotoEMCalPrimaryGeneratorAction();
  primaryGenerator->SetOptionConfiguration(options);
  runManager->SetUserAction(primaryGenerator);

  runManager->Initialize();

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  // G4UIExecutive* ui = new G4UIExecutive(argc, argv);
  if (options["mac"] == "vis.mac") {
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
  }
  if (options["mac"] != "") {
    G4String command = "/control/execute ";
    command += options["mac"];
    UImanager->ApplyCommand(command);
  } else if (options["nEvents"] != "") {
    G4String command = "/run/beamOn ";
    command += options["nEvents"];
    UImanager->ApplyCommand(command);
  }

  //  if (options["visualizationOnly"] != "true") {
  tf->cd();
  tr->Write(TObject::kOverwrite);
  tf->Close();
  //  }
  // delete visManager;

  return 0;
}

map<string, string> parseOptions(int argc, char** argv) {
  // Option parsing
  map<string, string> opt;
  vector<TString> argList;
  for (int iarg = 1; iarg < argc; iarg++) {
    string optString = string(argv[iarg]);
    istringstream iss(optString);
    string token;
    vector<string> tokens;
    while (iss >> token)
      tokens.push_back(token);
    for (auto tok : tokens) argList.push_back(tok.c_str());
  }

  for (int iarg = 0; iarg < argList.size(); iarg++) {
    if (argList[iarg].BeginsWith("-")) {
      argList[iarg].Remove(0, 1);
      if (argList[iarg].Contains("=")) {
        auto substr = argList[iarg].Tokenize("=");
        auto optionKey = ((TObjString*)substr->At(0))->GetString();
        auto optionValue = ((TObjString*)substr->At(1))->GetString();
        opt[optionKey.Data()] = optionValue.Data();
      } else {
        opt[argList[iarg].Data()] = argList[iarg + 1].Data();
        iarg++;
      }
    }
  }
  return opt;
}