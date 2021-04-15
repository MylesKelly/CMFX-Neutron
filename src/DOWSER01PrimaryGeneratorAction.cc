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
  // nAngle = 0.;
  G4double Emax = 6; // 1 MeV;
  G4double Emin = -0.39794; // 0.4 eV;
  G4double nEnergyT(0);
  fnArray = 1295;
  fAngleArray = 180;
  G4double fEnergyT = 0;
  pi = 3.14159265358979323846264338328;
  // Moon radius: 1737 km, LRO altitude: orbitH km
  orbitH = 50.0; // ground (km)
  moonR = 1737.0;
  RoverRH = (moonR/(moonR + orbitH));
  // vLRO = sqrt(G * moon_Mass /((orbitH+1737)*1000)); // LEO velocity in Y direction
  vLRO = sqrt(4902808346033.4/((orbitH+moonR)*1000)); //
  nEqvEnergy = 0.5*nMass*vLRO*vLRO/1.6022E-19; // eV
  // neutron potential difference
  // del U = G * n_Mass * moon_Mass * (1 / R_moon - 1 / (R_moon + H))
  nPotDiff = 51252.715*(1./1737000 - 1./(1737000 + orbitH*1000)); // potential difference: U_orbit - U_surf
  nMass = 1.6749E-27; // kg
  moonMass = 7.34767309E+22; // kg
  thetaCr = asin(RoverRH);
  // nMC2 = nMass * 299790000.0 * 299790000.0; // joule
  
  std::ifstream theData("LNeutronIntSpectrum.DAT", std::ios::in);
  G4cout << "Open Neutron Spectrum Data " << G4endl;
  
  if(!theData)
  {
    theData.close();
    G4cout << "LNeutronIntSpectrum.DAT does not exit! " << G4endl;
    G4double dElog10 = (Emax - Emin)/fnArray;
    for (G4int copyNo=0; copyNo<fnArray; copyNo++)
    {
      nEnergy[copyNo]=pow(10,Emin + copyNo * dElog10);
      intSpectrum[copyNo] = pow(nEnergy[copyNo],-0.89)
      *(pow(10,Emin + (copyNo+0.5) * dElog10)-pow(10,Emin + (copyNo-0.5) * dElog10));
      nEnergyT += intSpectrum[copyNo];
      intSpectrum[copyNo] = nEnergyT;
    }
    intSpectrum[0] = 0;
    for (G4int copyNo=0; copyNo<fnArray; copyNo++)
    {
      intSpectrum[copyNo] /=  nEnergyT;
    }
    
  } else
  {
    for (G4int copyNo=0; copyNo<fnArray; copyNo++)
    {
      if (!theData.eof())
      {
        theData >> nEnergy[copyNo] >> intSpectrum[copyNo]  ;
      }
    }
  }
  theData.close();
  
  // Lunar neutron angular distribution, LEND's View
  fEnergyT = 0;
  for (G4int copyNo = 0; copyNo < fAngleArray; copyNo++)
  {
    // nFluxHEE[copyNo] = cos(0.5*copyNo*pi/fAngleArray)*std::sqrt(cos(0.5*copyNo*pi/fAngleArray))*sin(0.5*copyNo*pi/fAngleArray);
    nFluxHEE[copyNo] = std::pow(cos(0.5*copyNo*pi/fAngleArray), 3./2.)*sin(0.5*copyNo*pi/fAngleArray);
    nThetaHEE[copyNo] = copyNo*90./fAngleArray;
    fEnergyT += nFluxHEE[copyNo];
    nFluxHEE[copyNo] = fEnergyT;
  }
  
  for (G4int copyNo = 0; copyNo < fAngleArray; copyNo++)
  {
    nFluxHEE[copyNo] /=  fEnergyT;
  }
  
  
  G4ParticleDefinition* particleDefinition
  = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
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

void DOWSER01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event
  G4double Eneutron(0), x0(0), y0(0), z0(0), px0(0), py0(0), pz0(0), px1(0), pz1(0), rotationAngle(0), theta(0), phi(0), x1(0), y1(0), 
  z1(0), Al_z(0), Xe_z(0), sourceRadius(0), HDPE_z(0), spread(0), spread_angle(0);
  G4int copyNo(0), iflag(0);
  G4double xx(0), yy(0), zz(0), pPhiH(0), pVerH(0), xDiff(0), yDiff(0), zDiff(0), U_torus(0), V_torus(0), C(0), A(0);
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
  
  //
  Al_z = 0.8 *mm;
  Xe_z = 7 *mm;
  HDPE_z = 40*mm;

  sourceRadius = 60*cm;  //+ Xe_z/2 + Al_z ;

  //Randomly generate rotation angle
  rotationAngle = pi * G4UniformRand() + pi/2;

  U_torus = 2 * pi * G4UniformRand();
  V_torus = 2 * pi * G4UniformRand();
  C = 1;
  A = 0.5;

  //Cartesian coordinate representation of a Toroidal surface following the wolfram website definition
  x1 = (C + A * cos(V_torus)) * cos(U_torus);
  y1 = (C + A * cos(V_torus)) * sin(U_torus);
  z1 = A * sin (V_torus);

  //Randonly generate spherical angle for momentum in specified range
  spread_angle = 6 * pi/180;
  G4double max_V = (cos(spread_angle) + 1)/2;

  G4double U = G4UniformRand();
  G4double V = G4UniformRand() * (1 - max_V);
  theta = 2*pi*U;
  phi = acos(2*(1 - V) - 1) ; //-spread_angle + 2*spread_angle*G4UniformRand();

  //generate cartesian coordinates for momentum
  px0 = 0;
  py0 = 0;
  pz0 = 1;

  xDiff = (cos(rotationAngle)) - x1;
  yDiff =  - y1;
  zDiff = (sin(rotationAngle)) - z1;

  fParticleGun->SetParticlePosition(G4ThreeVector(x1, y1, z1));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px0,py0,pz0));
  info->SetTheta0(rotationAngle*180.0/pi);
  info->SetTheta1(phi*180.0/pi);
  //analysisManager->FillH1(2, rotationAngle*180.0/pi);
  fParticleGun->GeneratePrimaryVertex(anEvent);
  anEvent->SetUserInformation(info);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


