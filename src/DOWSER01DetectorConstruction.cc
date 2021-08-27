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
/// \file DOWSER01DetectorConstruction.cc
/// \brief Implementation of the DOWSER01DetectorConstruction class

#include "DOWSER01DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Ellipsoid.hh"
#include "G4Para.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4SDManager.hh"
#include "G4SDChargedFilter.hh"
#include "G4SDParticleFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSEnergyDeposit3D.hh"
#include "G4PSTrackLength.hh"

#include "G4tgbRotationMatrix.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4GenericMessenger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <stdio.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01DetectorConstruction::DOWSER01DetectorConstruction()
: G4VUserDetectorConstruction(),
fCheckOverlaps(true)
{ ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01DetectorConstruction::~DOWSER01DetectorConstruction()
{ ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DOWSER01DetectorConstruction::Construct()
{
  // MATERIAL DEFINITIONS:
  // Setting up the NIST Manager object used for referencing the GEANT4 Material Databse:
  G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;

  //Specify material to serve as vacuum
  nistManager->FindOrBuildMaterial("G4_Galactic", fromIsotopes);
  G4Material* vacuum  = G4Material::GetMaterial("G4_Galactic");
  
  //Manually specifying Boron-10 material, allowing for precise control of denisty and atomic mass etc.
  G4int protons = 5, neutrons = 5, nucleons=protons+neutrons, isotopes, ncomponents;
  G4double atomicMass = 10.01294*g/mole;
  //Create Boron-10 Isotope:
  G4Isotope* isoB10 = new G4Isotope("isoB10", protons, nucleons, atomicMass);
  //Create Boron Element object:
  G4Element* elmB10 = new G4Element("elmB10", "B10", isotopes = 1);
  //Create Boron Material object
  G4Material* boron10 = new G4Material("boron10", 2.34*g/cm3, ncomponents=1, kStateSolid);

  //Specify that Boron Element consists of 100% Boron-10 Isotope
  elmB10->AddIsotope(isoB10,100*perCent);

  //Specify that Boron Material consists of 100% Boron Element
  boron10->AddElement(elmB10,100*perCent);
  
  //Define Xenon gas explicitly to allow for variable pressure, following same pattern as Boron above:
  G4Isotope* isoXe = new G4Isotope("isoXe131", 54, 131, 130.905*g/mole);
  G4Element* elmXe = new G4Element("elmXe", "XeGas", 1);
  elmXe->AddIsotope(isoXe, 100*perCent);

  //Note specification of temperature and pressure in the instaniation of the xenonGas material object:
  G4Material* xenonGas = new G4Material("xenonGas", 0.005894*g/cm3, ncomponents=1, kStateGas, 273*kelvin, 1*atmosphere);
  xenonGas->AddElement(elmXe, 100*perCent);

  //Load polyethylene material used as moderator from GEANT4 NIST database:
  nistManager->FindOrBuildMaterial("G4_POLYETHYLENE", fromIsotopes);
  G4Material* polyethylene = G4Material::GetMaterial("G4_POLYETHYLENE");

  // Print materials
  G4cout << G4endl << "The materials defined are: " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  G4cout << "************  end of material table  ************" << G4endl;
  //============================================================================
  //      Definitions of Solids, Logical Volumes, Physical Volumes
  //============================================================================
  
  //-------------
  // World Volume
  //-------------
  
  //Specify world volume:
  //Specify solid:
  G4ThreeVector worldSize = G4ThreeVector(150*cm, 150*cm, 150*cm);
  G4Box * solidWorld = new G4Box("soildWorld", worldSize.x()/2., worldSize.y()/2., worldSize.z()/2.);

  //Specify world logical volume:
  G4LogicalVolume * World = new G4LogicalVolume(solidWorld, vacuum, "World", 0, 0, 0);
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  G4VPhysicalVolume * worldPV
  = new G4PVPlacement(0,               // no rotation
                      G4ThreeVector(), // at (0,0,0)
                      World,      // its logical volume
                      "WorldPV",         // its name
                      0,               // its mother  volume
                      false,           // no boolean operations
                      0);              // copy number

  //DIMENSIONS FOR OBJECTS IN WORLD VOLUME:

  //Prism Parameters
  //The height of the isosceles triangle cross section:
  G4double prism_h = 175*mm;
  //1/2 of the base of the trangle cross section:
  G4double prism_l = (150/2)*mm;

  //HDPE Wedge:
  G4double HDPE_x1 = 2*prism_l;
  G4double HDPE_x2 = 0*mm;
  G4double HDPE_y1 = 100*mm;
  G4double HDPE_y2 = HDPE_y1;
  G4double HDPE_z = prism_h;

  //Detector Solid (Boron-10 Film):
  G4double Det_x = sqrt(pow(prism_l, 2) + pow(prism_h, 2));
  G4double Det_y = HDPE_y1;
  G4double Det_z = 0.001*mm;

  //SOLIDS:


  //G4VSolid* Attachment = new G4SubtractionSolid("Attachment", HDPESolid, HDPESubtraction, rotMat, G4ThreeVector(0, 0, HDPE_z/2));
  G4VSolid* HDPESolid = new G4Trd("HDPE", HDPE_x1/2, HDPE_x2/2, HDPE_y1/2, HDPE_y2/2, HDPE_z/2);

  G4VSolid* Detector = new G4Box("Detector", Det_x/2, Det_y/2, Det_z);

  //LOGICAL VOLUMES:
  G4LogicalVolume* HDPELogical = new G4LogicalVolume(HDPESolid, polyethylene, "HDPE");
  G4LogicalVolume* Detector1Logical = new G4LogicalVolume(Detector, boron10, "Detector1");
  G4LogicalVolume* Detector2Logical = new G4LogicalVolume(Detector, boron10, "Detector2");

  //Place volumes in the world:
  new G4PVPlacement(0, G4ThreeVector(0, 0, HDPE_z/2), HDPELogical, "HDPE", World, false, 0);

  //Compute rotation matricies to rotate Boron-10 films to align with the diagonal faces of the HDPE wedge
  G4RotationMatrix* rotMatPlus = new G4RotationMatrix();
  rotMatPlus->rotateY(-atan(prism_l/prism_h) - pi/2);

  G4RotationMatrix* rotMatMinus = new G4RotationMatrix();
  rotMatMinus->rotateY(atan(prism_l/prism_h) + pi/2);

  G4VPhysicalVolume* Detector1Physical = new G4PVPlacement(rotMatMinus, G4ThreeVector(prism_l/2 + (Det_x/prism_h)*(Det_z), 0, 0), Detector1Logical, "Detector1", HDPELogical, false, 0);
  G4VPhysicalVolume* Detector2Physical = new G4PVPlacement(rotMatPlus, G4ThreeVector(-(prism_l/2 + (Det_x/prism_h)*(Det_z)), 0, 0), Detector2Logical, "Detector2", HDPELogical, false, 0);

  // Visualization attributes
  G4VisAttributes* whiteBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes* greenBoxVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,0.0,0.25));
  G4VisAttributes* redBoxVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0,0.75));
  G4VisAttributes* blueBoxVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.25));
  G4VisAttributes* yellowBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.0,0.25));
  G4VisAttributes* lightblueVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,0.5,0.15));
  G4VisAttributes* grayVisAtt= new G4VisAttributes(G4Colour(0.3,0.3,0.3,0.1));

  HDPELogical->SetVisAttributes (whiteBoxVisAtt);
  Detector1Logical->SetVisAttributes (redBoxVisAtt);
  Detector2Logical->SetVisAttributes (blueBoxVisAtt);
  
  // Always return the physical World
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

