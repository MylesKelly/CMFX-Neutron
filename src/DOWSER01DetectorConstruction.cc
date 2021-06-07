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
  //=====================
  // Material Definitions
  //=====================
  //
  //-------- NIST Materials ----------------------------------------------------
  //  Material Information imported from NIST database.
  //
  // Lead material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;
  // nistManager->FindOrBuildMaterial("G4_Pb", fromIsotopes);
  nistManager->FindOrBuildMaterial("G4_Galactic", fromIsotopes);
  G4Material* vacuum  = G4Material::GetMaterial("G4_Galactic");
  
  // nistManager->FindOrBuildMaterial("G4_B", fromIsotopes);
  // G4Material* boron10 = G4Material::GetMaterial("G4_B");
  // density = 1.0 *g/cm3;
  // a = 10.01294*g/mole;
  // G4Material* boron10 = new G4Material("boron10", 5., 10.01294*g/mole, 1.0*g/cm3);
  // nistManager->FindOrBuildMaterial("boron10", fromIsotopes);
  
  G4int protons = 5, neutrons = 5, nucleons=protons+neutrons, isotopes, ncomponents;
  G4double atomicMass = 10.01294*g/mole;
  G4Isotope* isoB10 = new G4Isotope("isoB10", protons, nucleons, atomicMass);
  G4Element* elmB10 = new G4Element("elmB10", "B10", isotopes = 1);
  elmB10->AddIsotope(isoB10,100*perCent);
  G4Material* boron10 = new G4Material("boron10", 2.34*g/cm3, ncomponents=1, kStateSolid);
  boron10->AddElement(elmB10,100*perCent);
  
  //Define Xenon gas explicitly to allow for variable pressure:
  G4Isotope* isoXe = new G4Isotope("isoXe131", 54, 131, 130.905*g/mole);
  G4Element* elmXe = new G4Element("elmXe", "XeGas", 1);
  elmXe->AddIsotope(isoXe, 100*perCent);

  G4Material* xenonGas = new G4Material("xenonGas", 0.005894*g/cm3, ncomponents=1, kStateGas, 273*kelvin, 1*atmosphere);
  xenonGas->AddElement(elmXe, 100*perCent);

  nistManager->FindOrBuildMaterial("G4_POLYETHYLENE", fromIsotopes);
  G4Material* polyethylene = G4Material::GetMaterial("G4_POLYETHYLENE");

  G4Material* Silicon = nistManager->FindOrBuildMaterial("G4_Si", fromIsotopes);
  //at STP
  G4Material* Xenon = nistManager->FindOrBuildMaterial("G4_Xe", fromIsotopes);

  // Al for substract, Aluminum
  G4Material* Aluminum = nistManager->FindOrBuildMaterial("G4_Al", fromIsotopes);

  /*
   // Vacuum
   new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
   kStateGas, 2.73*kelvin, 3.e-18*pascal);
   */
  G4String symbol;
  G4double a, z, density;


  //Cadmium
  G4Material* cadmium = nistManager->FindOrBuildMaterial("G4_Cd", true);
  G4Element* Cd113 = new G4Element("Cadmium113",symbol="Cd113", z= 48., a =112.904*g/mole);
  G4Material* cadmium113 = new G4Material("Cadmium113", density=8.7*g/cm3, ncomponents=1, kStateSolid);
  cadmium113->AddElement(Cd113, 1);
  
  // BC545 Materials
  // Natural Boron-loaded Premium Plastic Scintillator
  // 5% Boronncomponents=

  G4Element* C = new G4Element("Carbon", symbol="C", z=6., a= 12.01*g/mole);
  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  G4Element* B10 = new G4Element("Boron10",symbol="B10", z= 5., a =10.01294*g/mole);
  G4Material* BC545 = new G4Material("BC545", density=1.026*g/cm3, ncomponents=3, kStateSolid);
  G4double fractionmass;
  BC545->AddElement(C, fractionmass = 0.9018);
  BC545->AddElement(H, fractionmass = 0.0887);
  BC545->AddElement(B10, fractionmass = 0.0095);
  
  nistManager->FindOrBuildMaterial("G4_BGO", fromIsotopes);
  G4Material* BGO  = G4Material::GetMaterial("G4_BGO");
  
  // Boron Cabride
  // B4C
  G4int natoms;
  G4Material* B4C = new G4Material("B4C", density = 2.52*g/cm3, ncomponents = 2, kStateSolid);
  B4C->AddElement(C, natoms = 1);
  B4C->AddElement(B10, natoms = 4);
  
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
  
  G4ThreeVector worldSize = G4ThreeVector(150*cm, 150*cm, 150*cm);
  G4Box * solidWorld
  = new G4Box("soildWorld", worldSize.x()/2., worldSize.y()/2., worldSize.z()/2.);
  G4LogicalVolume * World
  = new G4LogicalVolume(solidWorld, vacuum, "World", 0, 0, 0);
  
  //
  //  Must place the World Physical volume unrotated at (0,0,0).
  G4VPhysicalVolume * worldPV
  = new G4PVPlacement(0,               // no rotation
                      G4ThreeVector(), // at (0,0,0)
                      World,      // its logical volume
                      "WorldPV",         // its name
                      0,               // its mother  volume
                      false,           // no boolean operations
                      0);              // copy number


  //Xenon Gas Dimensions:
  G4double Xenon_x = 25.4 *mm;
  G4double Xenon_y = 25.4 *mm;
  G4double Xenon_z = 7 *mm;

  //Aluminum Substrate:
  G4double Al_x = 100 *mm;
  G4double Al_y = 100 *mm;
  G4double Al_z = 0.8 *mm;

  //Boron-10 Film:
  G4double Boron_x = 100 *mm;
  G4double Boron_y = 100 *mm;
  G4double Boron_z = 0.001 *mm;

  //SiPM:
  G4double SiPM_x = 3*mm; 
  G4double SiPM_y = 1*mm; 
  G4double SiPM_z = 3*mm;

  //Prism Parameters
  //The height of the isosceles triangle cross section:
  G4double prism_h = 175*mm;
  //1/2 of the base of the trangle cross section:
  G4double prism_l = (150/2)*mm;

  //HDPE Attachment:
  G4double HDPE_x1 = 2*prism_l;
  G4double HDPE_x2 = 0*mm;
  G4double HDPE_y1 = 100*mm;
  G4double HDPE_y2 = HDPE_y1;
  G4double HDPE_z = prism_h;

  //HDPE Subtraction Solid:
  G4double Sub_x = 16*mm;
  G4double Sub_y = Al_y + 2*mm;
  G4double Sub_z = 16*mm;

  //Detector Solid:
  G4double Det_x = sqrt(pow(prism_l, 2) + pow(prism_h, 2));
  G4double Det_y = HDPE_y1;
  G4double Det_z = 0.001*mm;

  //SOLIDS:
  G4VSolid* XenonSolid = new G4Box("XenonGas", Xenon_x/ 2, Xenon_y/2, Xenon_z/2);
  G4VSolid* AlSubstrateSolid = new G4Box("AlSubstrate", Al_x/ 2, Al_y/2, Al_z/2);
  G4VSolid* BoronFilmSolid = new G4Box("BoronFilm", Boron_x/ 2, Boron_y/2, Boron_z/2);
  G4VSolid* SiPMSolid = new G4Box("SiPM", SiPM_x/2, SiPM_y/2, SiPM_z/2);

  //G4VSolid* HDPESolid = new G4Box("HDPE", HDPE_x1/2, HDPE_y1/2, HDPE_z/2);
  G4VSolid* HDPESubtraction = new G4Box("HDPE_Sub", Sub_x/2, Sub_y/2, Sub_z/2);

  G4RotationMatrix* rotMatPlus = new G4RotationMatrix();
  rotMatPlus->rotateY(-atan(prism_l/prism_h) - pi/2);

  G4RotationMatrix* rotMatMinus = new G4RotationMatrix();
  rotMatMinus->rotateY(atan(prism_l/prism_h) + pi/2);

  //G4VSolid* Attachment = new G4SubtractionSolid("Attachment", HDPESolid, HDPESubtraction, rotMat, G4ThreeVector(0, 0, HDPE_z/2));
  G4VSolid* HDPESolid = new G4Trd("HDPE", HDPE_x1/2, HDPE_x2/2, HDPE_y1/2, HDPE_y2/2, HDPE_z/2);

  G4VSolid* Detector = new G4Box("Detector", Det_x/2, Det_y/2, Det_z);

  //LOGICAL VOLUMES:
  G4LogicalVolume* XenonLogical = new G4LogicalVolume(XenonSolid, Xenon, "XenonGas");
  G4LogicalVolume* AlSubstrateLogical = new G4LogicalVolume(AlSubstrateSolid, Aluminum, "AlSubstrate");
  G4LogicalVolume* BoronFilmLogical = new G4LogicalVolume(BoronFilmSolid, boron10, "BoronFilm");
  G4LogicalVolume* SiPMLogical = new G4LogicalVolume(SiPMSolid, Silicon, "SiPM");

  G4LogicalVolume* HDPELogical = new G4LogicalVolume(HDPESolid, polyethylene, "HDPE");

  G4LogicalVolume* Detector1Logical = new G4LogicalVolume(Detector, boron10, "Detector1");
  G4LogicalVolume* Detector2Logical = new G4LogicalVolume(Detector, boron10, "Detector2");

  G4int nSlices = 11;

  G4VSolid* detSlice = new G4Box("detSlice", Det_x/(nSlices), Det_y/2, Det_z);
  G4LogicalVolume* detSliceLogical = new G4LogicalVolume(detSlice, boron10, "detSlice");

  

  //Place volumes in the world:
  // new G4PVPlacement(0, G4ThreeVector(), XenonLogical, "XenonGas", World, false, 0);
  // new G4PVPlacement(0, G4ThreeVector(0, 0, Xenon_z/2 + Al_z/2 + Boron_z), AlSubstrateLogical, "AlSubstrate", World, false, 0);
  // new G4PVPlacement(0, G4ThreeVector(0, 0, -Xenon_z/2 - Al_z/2), AlSubstrateLogical, "AlSubstrate", World, false, 1);
  // new G4PVPlacement(0, G4ThreeVector(0, 0, Xenon_z/2 + Boron_z/2), BoronFilmLogical, "BoronFilm", World, false, 0);

  // new G4PVPlacement(0, G4ThreeVector(Xenon_x/8, Xenon_y/2 + SiPM_y/2, 0), SiPMLogical, "SiPM", XenonLogical, false, 0);
  // new G4PVPlacement(0, G4ThreeVector(Xenon_x/8, -Xenon_y/2 - SiPM_y/2, 0), SiPMLogical, "SiPM", XenonLogical, false, 1);

  // new G4PVPlacement(0, G4ThreeVector(Xenon_x/8 + Xenon_x/4, Xenon_y/2 + SiPM_y/2, 0), SiPMLogical, "SiPM", XenonLogical, false, 2);
  // new G4PVPlacement(0, G4ThreeVector(Xenon_x/8 + Xenon_x/4, -Xenon_y/2 - SiPM_y/2, 0), SiPMLogical, "SiPM", XenonLogical, false, 3);

  // new G4PVPlacement(0, G4ThreeVector(-Xenon_x/8, Xenon_y/2 + SiPM_y/2, 0), SiPMLogical, "SiPM", XenonLogical, false, 4);
  // new G4PVPlacement(0, G4ThreeVector(-Xenon_x/8, -Xenon_y/2 - SiPM_y/2, 0), SiPMLogical, "SiPM", XenonLogical, false, 5);

  // new G4PVPlacement(0, G4ThreeVector(-Xenon_x/8 - Xenon_x/4, Xenon_y/2 + SiPM_y/2, 0), SiPMLogical, "SiPM", XenonLogical, false, 6);
  // new G4PVPlacement(0, G4ThreeVector(-Xenon_x/8 - Xenon_x/4, -Xenon_y/2 - SiPM_y/2, 0), SiPMLogical, "SiPM", XenonLogical, false, 7);

  new G4PVPlacement(0, G4ThreeVector(0, 0, Xenon_z/2 + Al_z + Boron_z + HDPE_z/2), HDPELogical, "HDPE", World, false, 0);

  G4VPhysicalVolume* Detector1Physical = new G4PVPlacement(rotMatMinus, G4ThreeVector(prism_l/2 + (Det_x/prism_h)*(Det_z), 0, 0), Detector1Logical, "Detector1", HDPELogical, false, 0);
  G4VPhysicalVolume* Detector2Physical = new G4PVPlacement(rotMatPlus, G4ThreeVector(-(prism_l/2 + (Det_x/prism_h)*(Det_z)), 0, 0), Detector2Logical, "Detector2", HDPELogical, false, 0);

  //G4VPhysicalVolume* Detector1_rep = new G4PVReplica("Detector1_rep", detSliceLogical, Detector1Logical, kXAxis, nSlices, Det_x/(nSlices)); 
  //G4VPhysicalVolume* Detector3_rep = new G4PVReplica("Detector3_rep", detSliceLogical, Detector2Logical, kXAxis, nSlices, Det_x/(nSlices)); 

  //
  // Visualization attributes
  //
  // WorldLV->SetVisAttributes (G4VisAttributes::Invisible);
  G4VisAttributes* whiteBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  // G4VisAttributes* redBoxVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  G4VisAttributes* greenBoxVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,0.0,0.25));
  G4VisAttributes* redBoxVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0,0.75));
  G4VisAttributes* blueBoxVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.25));
  G4VisAttributes* yellowBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.0,0.25));
  G4VisAttributes* lightblueVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,0.5,0.15));
  G4VisAttributes* grayVisAtt= new G4VisAttributes(G4Colour(0.3,0.3,0.3,0.1));

  BoronFilmLogical->SetVisAttributes (redBoxVisAtt);
  AlSubstrateLogical->SetVisAttributes (blueBoxVisAtt);
  XenonLogical->SetVisAttributes (greenBoxVisAtt);
  SiPMLogical->SetVisAttributes (yellowBoxVisAtt);
  HDPELogical->SetVisAttributes (whiteBoxVisAtt);
  Detector1Logical->SetVisAttributes (redBoxVisAtt);
  Detector2Logical->SetVisAttributes (blueBoxVisAtt);
  // Cad1Logical->SetVisAttributes (blueBoxVisAtt);
  // Cad2Logical->SetVisAttributes (blueBoxVisAtt);
  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

