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
// $Id$
//
/// \file DOWSER01RunAction.cc
/// \brief Implementation of the DOWSER01RunAction class

#include "DOWSER01RunAction.hh"
#include "DOWSER01Analysis.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01RunAction::DOWSER01RunAction()
 : G4UserRunAction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01RunAction::~DOWSER01RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DOWSER01RunAction::BeginOfRunAction(const G4Run* run)
{ 
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;

  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Book histograms, ntuple
  //
  
  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in DOWSER01Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() 
         << " analysis manager" << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  
  // Open an output file
  //
  G4String fileName = "CMFX_Neutron_Radius_Test";
  analysisManager->OpenFile(fileName);
  analysisManager->SetFirstHistoId(1);

  
  // int filmBin = 250;
  // double filmMax = 2.5; // micron, -> +/- 1250 nm,
  int filmBin = 300;
  double filmMax = 30.0; // mm, -> +/- 1250 nm,

  int rangeBin = 100;
  double rangeMax = 70.0, rangeMin = -rangeMax; // mm
  
  // Creating histograms
  //
  // Lunar Neutrons:
  // 1 - energy
  // 2 - emission angle on lunar surface
  // 3 - entry angle at LRO altitude
  //analysisManager->CreateH1("nEnergy","HDPE Neutron Energy Spectrum; Log_{10} E(eV)", 120, -4., 8.);         // h1 1 incident N, E
  //analysisManager->CreateH1("nAngleG","HDPE Incident Neutron Angular Distribition; #theta0", 180, -90, 90);  // h1 2 incident N, theta
  //analysisManager->CreateH1("nDetEnergy1","Detected Neutron Spectrum; Log_{10} E(eV)", 120, -4., 8.);        // h1 3 detected N, E
  //analysisManager->CreateH1("nDetAngleG1","Detected Neutron Angular Distribition; #theta", 180, -90, 90);   // h1 4 detected N, theta
  //analysisManager->CreateH1("alphaEnergyB10","HDPE #alpha Spectrum; E(MeV)", 100,0,4.0);                    // h1 5 detected alpha, E
  //analysisManager->CreateH1("alphaFinPos1","HDPE #alpha Final Position Total, Z (mm)", rangeBin, rangeMin, rangeMax);// h1 6 alpha, final Z
  // analysisManager->CreateH1("li7nergyB10","HDPE Li7 Spectrum; E(MeV)", 100,0,4.0);                         // h1 7 detected Li, E
  // analysisManager->CreateH1("LiFinalPos1","HDPE Li7 Final Position Total, Z (mm)", rangeBin, rangeMin, rangeMax);  // h1 8 Li, final 
  // analysisManager->CreateH1("alphaEnergyXe","HDPE #alpha Spectrum Entering Xe; E(MeV)", 100,0,4.0);    // h1 9 detected alpha, E
  // analysisManager->CreateH1("li7EnergyXe","HDPE Li7 Spectrum Entering Xe; E(MeV)", 100,0,4.0);         // h1 10 detected Li, E
  // analysisManager->CreateH2("alphaFinalPos","HDPE #alpha Final Position in Xe(mm), X-Z", rangeBin, rangeMin, rangeMax, rangeBin, rangeMin, rangeMax);    // h2 1
  // analysisManager->CreateH2("LiFinalPos","HDPE Li7 Final Position in Xe (mm), X-Z", rangeBin, rangeMin, rangeMax, rangeBin, rangeMin, rangeMax);         // h2 2
  // analysisManager->CreateH2("thetaEnergy","Theta-Energy Log_{10} E(eV)", 180, -90, 90, 120, -4., 8.);    // h2 3
  analysisManager->CreateH1("alphaEnergy", "Alpha Particle Energy Spectrum in MeV", 10, 1., 2.);
  analysisManager->CreateH1("nEnergyB10", "Neutron Energy incident on Boron-10 Thin Film", 50, 0.000000001, 2500000, "none", "none", "log");
  analysisManager->CreateH1("b10rotation", "Boron-10 Count vs. Rotation Angle", 25, 0, 180);
  analysisManager->CreateH1("b10phi", "Boron-10 Count vs. Phi", 50, -5, 5);
  analysisManager->CreateH1("D1b10", "Boron-10 Count vs. Rotation Angle in Detector 1", 30, 0, 180);
  analysisManager->CreateH1("D2b10", "Boron-10 Count vs. Rotation Angle in Detector 2", 30, 0, 180);
  analysisManager->CreateH1("D1Li_x", "Boron-10 Count vs. Detector Depth", 50, 0, 200);

  analysisManager->CreateH1("D1b10_1", "Boron-10 Count vs. Rotation Angle in Detector 1, Segment 1", 30, 0, 180);
  analysisManager->CreateH1("D2b10_1", "Boron-10 Count vs. Rotation Angle in Detector 2, Segment 1", 30, 0, 180);
  analysisManager->CreateH1("D1b10_2", "Boron-10 Count vs. Rotation Angle in Detector 1, Segment 2", 30, 0, 180);
  analysisManager->CreateH1("D2b10_2", "Boron-10 Count vs. Rotation Angle in Detector 2, Segment 2", 30, 0, 180);
  analysisManager->CreateH1("D1b10_3", "Boron-10 Count vs. Rotation Angle in Detector 1, Segment 3", 30, 0, 180);
  analysisManager->CreateH1("D2b10_3", "Boron-10 Count vs. Rotation Angle in Detector 2, Segment 3", 30, 0, 180);
  analysisManager->CreateH1("D1b10_4", "Boron-10 Count vs. Rotation Angle in Detector 1, Segment 4", 30, 0, 180);
  analysisManager->CreateH1("D2b10_4", "Boron-10 Count vs. Rotation Angle in Detector 2, Segment 4", 30, 0, 180);

  analysisManager->CreateH1("nRadius", "Neutron Production vs. Radius", 150, 0, 500);

  analysisManager->CreateH2("nEnergyB10rotationD1", "Incoming Neutron Energy Spectrum on Boron-10 as a function of rotationAngle, Detector 1", 25, 0, 180, 100, 0.001, 1, "none", "none", "none", "none", "linear","log");
  analysisManager->CreateH2("nEnergyB10rotationD2", "Incoming Neutron Energy Spectrum on Boron-10 as a function of rotationAngle, Detector 2", 25, 0, 180, 100, 0.001, 1, "none", "none", "none", "none", "linear","log");
  analysisManager->CreateH2("Li_xB10rotation", "Boron 10 Capture X position vs. Rotation Angle", 25, 0, 180, 50, 0, 200);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DOWSER01RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int nofEvents = aRun->GetNumberOfEvent();
  if ( nofEvents == 0 ) return;
  
  // print histogram statistics
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();


  // save histograms 
  //
  analysisManager->Write();
  analysisManager->CloseFile();
  
  // complete cleanup
  //
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
