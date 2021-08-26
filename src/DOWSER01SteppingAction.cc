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
/// \file DOWSER01SteppingAction.cc
/// \brief Implementation of the DOWSER01SteppingAction class

#include "DOWSER01SteppingAction.hh"
#include "DOWSER01DetectorConstruction.hh"
#include "DOWSER01Analysis.hh"  // 2012-0403, J.J. Su
#include "DOWSER01EventInformation.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4CsvAnalysisManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4GenericMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01SteppingAction* DOWSER01SteppingAction::fgInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01SteppingAction* DOWSER01SteppingAction::Instance()
{
    // Static acces function via G4RunManager
    
    return fgInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01SteppingAction::DOWSER01SteppingAction()
: G4UserSteppingAction(), kEnergy(0), nEnergy(0), fParticleName(), numOfCapture(0)
{
    fgInstance = this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01SteppingAction::~DOWSER01SteppingAction()
{
    fgInstance = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DOWSER01SteppingAction::UserSteppingAction(const G4Step* step)
{
  // Get name of logical volume this step is in
  G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4String nameLogicVolume = volume->GetName();

  //Name of particle execxuting this step:
  fParticleName = step->GetTrack()->GetDefinition()->GetParticleName();

  //Get the number of the current step:
  G4int fStepNo = step->GetTrack()->GetCurrentStepNumber();
    
  //Get reference to analysisManger object so data can be stored in histograms:
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  if (fStepNo == 1)
  {
    G4EventManager* evtMan = G4EventManager::GetEventManager();
    DOWSER01EventInformation* info = (DOWSER01EventInformation*)evtMan->GetConstCurrentEvent()->GetUserInformation();
    if (info==0){G4cout<<"Warning: No Event information recorded!"<<G4endl;}
    
    else {
      theta0 = info->GetTheta0();
      nEnergy = info->GetEnergyN();
    }

    //kEnergy = step->GetPreStepPoint()->GetKineticEnergy()/CLHEP::eV;

    if (fParticleName == "Li7")
    {
      //Incremet Boron-10 capture counter when Li-7 ion is produced
      numOfCapture += 1;
      G4cout << "New B-10 capture event, total is: " << numOfCapture << G4endl;
    }
  } 

  //Get particle energy as it enters the Boron-10 thin film:
  if (nameLogicVolume == "Detector" && fParticleName == "neutron")
  {
    G4double stepEnergy = step->GetPostStepPoint()->GetKineticEnergy()/CLHEP::eV;

    //If Energy is not zero, store it in a histogram:
    if (stepEnergy != 0)
    {
        analysisManager->FillH1(2, stepEnergy);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DOWSER01SteppingAction::Reset()
{
    G4cout << "Stepping Action Reset Called" << G4endl;
}

