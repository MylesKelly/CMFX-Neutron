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
  //Log run number to the G4 console
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
  
  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in DOWSER01Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() 
         << " analysis manager" << G4endl;

  // Open an output file
  G4String fileName = "DOWSER_CMFX_DD_150x200mm_Segmented_Isotropic_100mil";
  analysisManager->OpenFile(fileName);
  analysisManager->SetFirstHistoId(1);

  // Creating histograms
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
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();


  // save histograms 
  analysisManager->Write();
  analysisManager->CloseFile();
  
  // complete cleanup
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
