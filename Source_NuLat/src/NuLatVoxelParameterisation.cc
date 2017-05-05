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
// $Id: NuLatVoxelParameterisation.cc 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file NuLatVoxelParameterisation.cc
/// \brief Implementation of the NuLatVoxelParameterisation class

#include "NuLatVoxelParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "stdlib.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NuLatVoxelParameterisation::NuLatVoxelParameterisation(
    G4int    nOfVoxelsInX,   G4int    nOfVoxelsInY,   G4int    nOfVoxelsInZ,
    G4double voxelXDim,      G4double voxelYDim,      G4double voxelZDim,
    G4double voxSpacingXDim, G4double voxSpacingYDim, G4double voxSpacingZDim)
          : G4VPVParameterisation()
{
  G4double totalNumberOfVoxels = nOfVoxelsInX*nOfVoxelsInY*nOfVoxelsInZ;
  xVoxel = (G4double*) malloc(totalNumberOfVoxels*sizeof(G4double));
  yVoxel = (G4double*) malloc(totalNumberOfVoxels*sizeof(G4double));
  zVoxel = (G4double*) malloc(totalNumberOfVoxels*sizeof(G4double));

  for (G4int copyNo=0;copyNo<totalNumberOfVoxels;copyNo++)
  {
    G4int xCellIndex = 1+(copyNo % (nOfVoxelsInX*nOfVoxelsInY))/ nOfVoxelsInY;
    G4int yCellIndex = 1+copyNo % nOfVoxelsInY;
    G4int zCellIndex = 1+copyNo / (nOfVoxelsInX*nOfVoxelsInY);

    G4double voxSpacingX = voxelXDim+voxSpacingXDim;
	  G4double voxSpacingY = voxelYDim+voxSpacingYDim;
	  G4double voxSpacingZ = voxelZDim+voxSpacingZDim;

    if(nOfVoxelsInX%2 == 0)
      xVoxel[copyNo] = (((xCellIndex) - (nOfVoxelsInX/2+0.5)) * voxSpacingX)*2;
    else
     xVoxel[copyNo] = ((xCellIndex-1 - nOfVoxelsInX/2) * voxSpacingX)*2;

    if(nOfVoxelsInY%2 == 0)
      yVoxel[copyNo] = (((yCellIndex) - (nOfVoxelsInY/2+0.5)) * voxSpacingY)*2;
    else
      yVoxel[copyNo] = ((yCellIndex-1 - nOfVoxelsInY/2) * voxSpacingY)*2; //Note nOfVoxelsInY/2 truncates the division
      
    if(nOfVoxelsInZ%2 == 0)
      zVoxel[copyNo] = (((zCellIndex) - (nOfVoxelsInZ/2+0.5)) * voxSpacingZ)*2;
    else
      zVoxel[copyNo] = ((zCellIndex-1 - nOfVoxelsInZ/2) * voxSpacingZ)*2; //Note nOfVoxelsInZ/2 truncates the division
    
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NuLatVoxelParameterisation::~NuLatVoxelParameterisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NuLatVoxelParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const
{
  physVol->SetTranslation(G4ThreeVector(xVoxel[copyNo],yVoxel[copyNo],zVoxel[copyNo]));
  physVol->SetCopyNo(copyNo);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
