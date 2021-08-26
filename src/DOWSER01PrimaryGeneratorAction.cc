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
/// \file DOWSER01PrimaryGeneratorAction.cc
/// \brief Implementation of the DOWSER01PrimaryGeneratorAction class

#include "DOWSER01PrimaryGeneratorAction.hh"
#include "DOWSER01Analysis.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "DOWSER01EventInformation.hh" // Event Info

// #include "G4CsvAnalysisManager.hh"  // 2013-02-04

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01PrimaryGeneratorAction::DOWSER01PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
fParticleGun(0)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  pi = 3.14159265358979323846264338328;

  G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
  fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(0.001032720*eV);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01PrimaryGeneratorAction::~DOWSER01PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double RadialPDF(G4double x)
{
  return exp(-20*pow(x - 0.5, 2));
  //return -pow(2 * x - 1, 2) + 1;
}

G4double AxialPDF(G4double x)
{
  return exp(-10*pow(x - 0.5, 2));
}

G4double GenerateNumFromPDF(G4double (*PDF)(G4double), G4double weight = 1)
{
  while(true)
  { 
    G4double num = G4UniformRand();
    G4double prob = PDF(num) * weight; 
    if (G4UniformRand() <= prob)
    {
      return num;
    }
  }
}

void DOWSER01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event
  G4double Eneutron(0), x0(0), y0(0), z0(0), px0(0), py0(0), pz0(0), px1(0), pz1(0), rotationAngle(0), theta(0), phi(0), x1(0), y1(0), 
  z1(0), Al_z(0), Xe_z(0), sourceRadius(0), HDPE_z(0), spread(0), spread_angle(0);
  G4int copyNo(0), iflag(0);
  G4double xx(0), yy(0), zz(0), pPhiH(0), pVerH(0), xDiff(0), yDiff(0), zDiff(0), R(0), Phi(0), Z(0), maxRadius(0), minRadius(0), radialPick(0);
  // G4double xx(0), yy(0), zz(0), pPhiH(0), pVerH(0);
  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore
  //
  G4double worldZHalfLength = 0;
  G4LogicalVolume* worlLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
  G4Box* worldBox = 0;
  if ( worlLV) worldBox = dynamic_cast< G4Box*>(worlLV->GetSolid());
  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength();
  }
  else  {
    G4cerr << "World volume of box not found." << G4endl;
    G4cerr << "Perhaps you have changed geometry." << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
  }
  
  // begining of particle iteration - energy, position, angle assessment
  
  // NEW - SET EVENT INFORMATION
  // Required files:
  // DOWSER01EventInformation.hh
  
  DOWSER01EventInformation* info = new DOWSER01EventInformation();
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  // begining of energy assigment
  
  //Set neutron energy
  Eneutron = 2.45;
  
  fParticleGun->SetParticleEnergy(Eneutron*MeV);
  analysisManager->FillH1(1, std::log10(Eneutron));
  info->SetEnergyN(Eneutron);
  
  maxRadius = 500;
  minRadius = 10;

  //Sqrt present to correct for point clustering in disk point picking
  //for more detail: https://mathworld.wolfram.com/DiskPointPicking.html
  radialPick = GenerateNumFromPDF(RadialPDF);
  R = maxRadius*sqrt(radialPick + pow((minRadius/maxRadius), 2));
  Phi = 2 * pi * G4UniformRand();
  Z = 300 * (2*GenerateNumFromPDF(AxialPDF, RadialPDF(radialPick))-1);

  analysisManager->FillH1(16, R);
  analysisManager->FillH2(4, R, Z);

  //Cylindrical coordinate representation of a Cylindrical surface following the wolfram website definition
  x1 = R * sin(Phi);
  y1 = R * cos(Phi);
  z1 = Z;

  //generate cartesian coordinates for momentum
  px0 = 1 - 2*G4UniformRand();
  py0 = 1 - 2*G4UniformRand(); 
  pz0 = 1 - 2*G4UniformRand(); 

  fParticleGun->SetParticlePosition(G4ThreeVector(x1, y1, z1));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px0,py0,pz0));
  info->SetTheta0(rotationAngle*180.0/pi);
  info->SetTheta1(phi*180.0/pi);
  //analysisManager->FillH1(2, rotationAngle*180.0/pi);
  fParticleGun->GeneratePrimaryVertex(anEvent);
  anEvent->SetUserInformation(info);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


