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

  //Load Stainless Steel for use in vacuum vessel:
  G4Material* Steel = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL", fromIsotopes);


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
  G4RotationMatrix* rot = new G4RotationMatrix();
  rot->rotateZ(0.*deg);
  
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

  //HDPE Attachment:
  G4double HDPE_x = 75*mm;
  G4double HDPE_y = 124*mm;
  G4double HDPE_z = 124*mm;

  //Detector Solid:
  G4double Det_x = 0.001*mm;
  G4double Det_y = HDPE_y;
  G4double Det_z = HDPE_z;

  //Boron-10 Film:
  G4double Boron_x = 0.001 *mm;
  G4double Boron_y = HDPE_y;
  G4double Boron_z = HDPE_z;

  //Vaccum Vessel
  G4double InnerRadius = (100-1.3)/2*cm;
  G4double OuterRadius = 100/2*cm;
  G4double length = 100*cm;

  //Central Conductor
  G4double conductorRadius = 5*cm;

  //HDPE Ring
  G4double HDPE_InnerRadius = OuterRadius + 5*cm;
  G4double HDPE_OuterRadius = HDPE_InnerRadius + 5*cm;
  G4double HDPE_length = length/2;

  //Detector Tube:
  G4double Det_InnerRadius = HDPE_OuterRadius;
  G4double Det_OuterRadius = Det_InnerRadius + 0.001*mm;
  G4double Det_length = HDPE_length;

  //SOLIDS:
  G4VSolid* BoronFilmSolid = new G4Box("BoronFilm", Boron_x/ 2, Boron_y/2, Boron_z/2);
  G4VSolid* VVesselSolid = new G4Tubs("VVessel", InnerRadius, OuterRadius, length/2, 0, 2*pi);
  G4VSolid* CConductorSolid = new G4Tubs("CConductor", 0, conductorRadius, length/2, 0, 2*pi);
  G4VSolid* HDPESolid = new G4Box("HDPESolid", HDPE_x/2, HDPE_y/2, HDPE_z/2);

  //G4VSolid* HDPERing = new G4Tubs("HDPERing", HDPE_InnerRadius, HDPE_OuterRadius, HDPE_length/2, 0, pi/2);
  //G4VSolid* RingDetector = new G4Tubs("Detector", Det_InnerRadius, Det_OuterRadius, Det_length/2, 0, pi/2);

  //LOGICAL VOLUMES:
  G4LogicalVolume* BoronFilmLogical = new G4LogicalVolume(BoronFilmSolid, boron10, "BoronFilm");
  G4LogicalVolume* HDPELogical = new G4LogicalVolume(HDPESolid, polyethylene, "HDPE");

  G4LogicalVolume* VVesselLogical = new G4LogicalVolume(VVesselSolid, Steel, "VVessel");  
  G4LogicalVolume* CConductorLogical = new G4LogicalVolume(CConductorSolid, Steel, "CConductor");  

  //G4LogicalVolume* RingDetectorLogical = new G4LogicalVolume(RingDetector, boron10, "RingDetector");
  //G4LogicalVolume* HDPERingLogical = new G4LogicalVolume(HDPERing, polyethylene, "HDPERing");


  //PHYSICAL VOLUMES:

  G4VPhysicalVolume* HDPEPhysical = new G4PVPlacement(0, G4ThreeVector(0.505*m + HDPE_x/2, 0, 0), HDPELogical, "HDPE", World, false, 0);
  G4VPhysicalVolume* Detector1Physical = new G4PVPlacement(0, G4ThreeVector(HDPE_x/2 + Boron_x/2, 0, 0), BoronFilmLogical, "Film 1", HDPELogical, false, 0);
  G4VPhysicalVolume* VVesselPhysical = new G4PVPlacement(0, G4ThreeVector(), VVesselLogical, "VVessel", World, false, 0);

  //G4VPhysicalVolume* RingDetectorPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), RingDetectorLogical, "Detector", World, false, 0)
  //G4VPhysicalVolume* HDPERingPhysical = new G4PVPlacement(0, G4ThreeVector(), HDPERingLogical, "HDPERing", World, false, 0);


  // VISUALIZATION ATTRIBUTES:
  G4VisAttributes* whiteBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes* greenBoxVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,0.0,0.25));
  G4VisAttributes* redBoxVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0,0.75));
  G4VisAttributes* blueBoxVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.25));
  G4VisAttributes* yellowBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.0,0.25));
  G4VisAttributes* lightblueVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,0.5,0.15));
  G4VisAttributes* grayVisAtt= new G4VisAttributes(G4Colour(0.3,0.3,0.3,0.1));

  BoronFilmLogical->SetVisAttributes (redBoxVisAtt);
  HDPELogical->SetVisAttributes (whiteBoxVisAtt);
  VVesselLogical->SetVisAttributes (whiteBoxVisAtt);
  //RingDetectorLogical->SetVisAttributes (redBoxVisAtt);
  //HDPERingLogical->SetVisAttributes (whiteBoxVisAtt);

  // Always return the World physical volume:
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

