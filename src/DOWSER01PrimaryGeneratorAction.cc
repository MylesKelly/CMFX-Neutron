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
  //Instantiate particle gun object:
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);
  
  //Useful value for later:
  pi = 3.14159265358979323846264338328;

  //Specify type of particle to be launched from the gun:
  G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
  fParticleGun->SetParticleDefinition(particleDefinition);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01PrimaryGeneratorAction::~DOWSER01PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DOWSER01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event
  G4double Eneutron(0), x0(0), y0(0), z0(0), px0(0), py0(0), pz0(0), px1(0), pz1(0), rotationAngle(0), theta(0), phi(0), x1(0), y1(0), 
  z1(0), Al_z(0), Xe_z(0), sourceRadius(0), HDPE_z(0), spread(0), spread_angle(0);
  G4int copyNo(0), iflag(0);
  G4double xx(0), yy(0), zz(0), pPhiH(0), pVerH(0), xDiff(0), yDiff(0), zDiff(0);

  // begining of particle iteration - energy, position, angle assessment
  
  // NEW - SET EVENT INFORMATION
  // Required files:
  // DOWSER01EventInformation.hh
  DOWSER01EventInformation* info = new DOWSER01EventInformation();
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //Set neutron energy
  Eneutron = 2.45;
  
  fParticleGun->SetParticleEnergy(Eneutron*MeV);
  analysisManager->FillH1(1, std::log10(Eneutron));
  info->SetEnergyN(Eneutron);
  
  //Distance from the origin of the neutron source
  sourceRadius = 60*cm;  

  //generate random rotation angle from uniform distribution
  rotationAngle = pi * G4UniformRand() + pi/2;

  //calculate cartesian coordinates for neutron origin
  x1 = sourceRadius * sin(rotationAngle); //-1.25 + 2.5 * G4UniformRand();
  y1 = 0; //-1.25 + 2.5 * G4UniformRand();
  z1 = sourceRadius * cos(rotationAngle) + Xe_z/2 + Al_z + HDPE_z/2;

  //Randonly generate spherical angle for momentum in specified range
  spread_angle = 6 * pi/180;
  G4double max_V = (cos(spread_angle) + 1)/2;

  //This business of U and V variables being first computed and then used to generate spherical angles
  //comes from spherical point picking, details described here: https://mathworld.wolfram.com/SpherePointPicking.html
  G4double U = G4UniformRand();
  G4double V = G4UniformRand() * (1 - max_V);
  theta = 2*pi*U;
  phi = acos(2*(1 - V) - 1) ; //-spread_angle + 2*spread_angle*G4UniformRand();

  //compute cartesian coordinates for momentum from spherical coordinates generated above
  px0 = sin(phi)*cos(theta);
  py0 = sin(phi)*sin(theta);
  pz0 = -cos(phi);

  //Rotate momentum vector about y axis through angle rotationAngle, so zero polar angle always points to the origin
  px1 = px0*cos(-rotationAngle) - pz0*sin(-rotationAngle);
  pz1 = px0*sin(-rotationAngle) + pz0*cos(-rotationAngle);

  fParticleGun->SetParticlePosition(G4ThreeVector(x1, y1, z1));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px1,py0,pz1));

  //Set parameters of event info object so the data can be accessed elsewhere, namely stepping action
  info->SetTheta0(rotationAngle*180.0/pi);
  info->SetTheta1(phi*180.0/pi);
  //analysisManager->FillH1(2, rotationAngle*180.0/pi);

  fParticleGun->GeneratePrimaryVertex(anEvent);
  anEvent->SetUserInformation(info);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


