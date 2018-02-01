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
// $Id: NuLatLightGuideParameterisation.cc 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file NuLatLightGuideParameterisation.cc
/// \brief Implementation of the NuLatLightGuideParameterisation class

#include "NuLatLightGuideParameterisation.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NuLatLightGuideParameterisation::NuLatLightGuideParameterisation(
    G4int    nOfVoxelsInX,   G4int    nOfVoxelsInY,   G4int    nOfVoxelsInZ,
    G4double voxelXDim,      G4double voxelYDim,      G4double voxelZDim,
    G4double voxSpacingXDim, G4double voxSpacingYDim, G4double voxSpacingZDim,
    G4int xselector, G4int yselector, G4int zselector)
          : G4VPVParameterisation()
{
  G4double totalNumberOfLightGuides = nOfVoxelsInX*nOfVoxelsInY*2
                               + nOfVoxelsInY*nOfVoxelsInZ*2
                               + nOfVoxelsInZ*nOfVoxelsInX*2;


  G4double voxSpacingX = voxelXDim+voxSpacingXDim;
  G4double voxSpacingY = voxelYDim+voxSpacingYDim;
  G4double voxSpacingZ = voxelZDim+voxSpacingZDim;
	  
  if(xselector == 0 && yselector == 0 && zselector == 1)
  {
    for(G4int copyNo=0;copyNo<totalNumberOfLightGuides;copyNo++)
    {
      if(copyNo<nOfVoxelsInX*nOfVoxelsInY)//nOfVoxelsInX*nOfVoxelsInY)
      {

        G4int xCellIndex = 1+(copyNo % (nOfVoxelsInX*nOfVoxelsInY))/ nOfVoxelsInY;
        G4int yCellIndex = 1+copyNo % nOfVoxelsInY;

        zLightGuide[copyNo] = 0;//(((nOfVoxelsInZ-1) - nOfVoxelsInZ/2) * voxSpacingZ)*2 + (1.378*2.54*cm+0.5*cm)/2;
        if(nOfVoxelsInX%2 == 1)
          xLightGuide[copyNo] = (((xCellIndex-1) - nOfVoxelsInX/2) * voxSpacingX)*2;
        else
          xLightGuide[copyNo] = (((xCellIndex) - (nOfVoxelsInX/2+0.5)) * voxSpacingX)*2;

        if(nOfVoxelsInY%2 == 1)
          yLightGuide[copyNo] = (((yCellIndex-1) - nOfVoxelsInY/2) * voxSpacingY)*2;
        else
          yLightGuide[copyNo] = (((yCellIndex) - (nOfVoxelsInY/2+0.5)) * voxSpacingY)*2;
          
        rotm[copyNo]  = new G4RotationMatrix(G4ThreeVector(1,0,0), 0);
      }
 
    }
  }
  else if(xselector == 0 && yselector == 0 && zselector == -1)
  {
    for(G4int copyNo=0;copyNo<totalNumberOfLightGuides;copyNo++)
    {
      if(copyNo<nOfVoxelsInX*nOfVoxelsInY)//nOfVoxelsInX*nOfVoxelsInY)
      {

        G4int xCellIndex = 1+(copyNo % (nOfVoxelsInX*nOfVoxelsInY))/ nOfVoxelsInY;
        G4int yCellIndex = 1+copyNo % nOfVoxelsInY;

        zLightGuide[copyNo] = 0;//(((nOfVoxelsInZ-1) - nOfVoxelsInZ/2) * voxSpacingZ)*2 + (1.378*2.54*cm+0.5*cm)/2;
        if(nOfVoxelsInX%2 == 1)
          xLightGuide[copyNo] = (((xCellIndex-1) - nOfVoxelsInX/2) * voxSpacingX)*2;
        else
          xLightGuide[copyNo] = (((xCellIndex) - (nOfVoxelsInX/2+0.5)) * voxSpacingX)*2;

        if(nOfVoxelsInY%2 == 1)
          yLightGuide[copyNo] = (((yCellIndex-1) - nOfVoxelsInY/2) * voxSpacingY)*2;
        else
          yLightGuide[copyNo] = (((yCellIndex) - (nOfVoxelsInY/2+0.5)) * voxSpacingY)*2;
        rotm[copyNo]  = new G4RotationMatrix(G4ThreeVector(1,0,0), 180*deg);
      }
 
    }
  }
  else if(xselector == 0 && yselector == 1 && zselector == 0)
  {
    for(G4int copyNo=0;copyNo<totalNumberOfLightGuides;copyNo++)
    {
      if(copyNo<nOfVoxelsInX*nOfVoxelsInZ)
      {

        G4int xCellIndex = 1+(copyNo % (nOfVoxelsInX*nOfVoxelsInZ))/ nOfVoxelsInZ;
        G4int zCellIndex = 1+copyNo % nOfVoxelsInZ;

        yLightGuide[copyNo] = 0;
        if(nOfVoxelsInX%2 == 1)
          xLightGuide[copyNo] = (((xCellIndex-1) - nOfVoxelsInX/2) * voxSpacingX)*2;
        else
          xLightGuide[copyNo] = (((xCellIndex) - (nOfVoxelsInX/2+0.5)) * voxSpacingX)*2;

        if(nOfVoxelsInZ%2 == 1)
          zLightGuide[copyNo] = (((zCellIndex-1) - nOfVoxelsInZ/2) * voxSpacingZ)*2;
        else
          zLightGuide[copyNo] = (((zCellIndex) - (nOfVoxelsInZ/2+0.5)) * voxSpacingZ)*2;
        rotm[copyNo]  = new G4RotationMatrix(G4ThreeVector(1,0,0), 90*deg);
      }
 
    }
  }
  else if(xselector == 0 && yselector == -1 && zselector == 0)
  {
    for(G4int copyNo=0;copyNo<totalNumberOfLightGuides;copyNo++)
    {
      if(copyNo<nOfVoxelsInX*nOfVoxelsInZ)
      {

        G4int xCellIndex = 1+(copyNo % (nOfVoxelsInX*nOfVoxelsInZ))/ nOfVoxelsInZ;
        G4int zCellIndex = 1+copyNo % nOfVoxelsInZ;

        yLightGuide[copyNo] = 0;
        if(nOfVoxelsInX%2 == 1)
          xLightGuide[copyNo] = (((xCellIndex-1) - nOfVoxelsInX/2) * voxSpacingX)*2;
        else
          xLightGuide[copyNo] = (((xCellIndex) - (nOfVoxelsInX/2+0.5)) * voxSpacingX)*2;

        if(nOfVoxelsInZ%2 == 1)
          zLightGuide[copyNo] = (((zCellIndex-1) - nOfVoxelsInZ/2) * voxSpacingZ)*2;
        else
          zLightGuide[copyNo] = (((zCellIndex) - (nOfVoxelsInZ/2+0.5)) * voxSpacingZ)*2;
        rotm[copyNo]  = new G4RotationMatrix(G4ThreeVector(1,0,0), 270*deg);
      }
 
    }
  }
  else if(xselector == 1 && yselector == 0 && zselector == 0)
  {
    for(G4int copyNo=0;copyNo<totalNumberOfLightGuides;copyNo++)
    {
      if(copyNo<nOfVoxelsInZ*nOfVoxelsInY)
      {

        G4int zCellIndex = 1+(copyNo % (nOfVoxelsInZ*nOfVoxelsInY))/ nOfVoxelsInY;
        G4int yCellIndex = 1+copyNo % nOfVoxelsInY;

        xLightGuide[copyNo] = 0;
        if(nOfVoxelsInX%2 == 1)
          yLightGuide[copyNo] = (((yCellIndex-1) - nOfVoxelsInY/2) * voxSpacingY)*2;
        else
          yLightGuide[copyNo] = (((yCellIndex) - (nOfVoxelsInY/2+0.5)) * voxSpacingY)*2;

        if(nOfVoxelsInZ%2 == 1)
          zLightGuide[copyNo] = (((zCellIndex-1) - nOfVoxelsInZ/2) * voxSpacingZ)*2;
        else
          zLightGuide[copyNo] = (((zCellIndex) - (nOfVoxelsInZ/2+0.5)) * voxSpacingZ)*2;
        rotm[copyNo]  = new G4RotationMatrix(G4ThreeVector(0,1,0), 270*deg);
      }
 
    }
  }
  else if(xselector == -1 && yselector == 0 && zselector == 0)
  {
    for(G4int copyNo=0;copyNo<totalNumberOfLightGuides;copyNo++)
    {
      if(copyNo<nOfVoxelsInZ*nOfVoxelsInY)
      {

        G4int zCellIndex = 1+(copyNo % (nOfVoxelsInZ*nOfVoxelsInY))/ nOfVoxelsInY;
        G4int yCellIndex = 1+copyNo % nOfVoxelsInY;

        xLightGuide[copyNo] = 0;
        if(nOfVoxelsInX%2 == 1)
          yLightGuide[copyNo] = (((yCellIndex-1) - nOfVoxelsInY/2) * voxSpacingY)*2;
        else
          yLightGuide[copyNo] = (((yCellIndex) - (nOfVoxelsInY/2+0.5)) * voxSpacingY)*2;

        if(nOfVoxelsInZ%2 == 1)
          zLightGuide[copyNo] = (((zCellIndex-1) - nOfVoxelsInZ/2) * voxSpacingZ)*2;
        else
          zLightGuide[copyNo] = (((zCellIndex) - (nOfVoxelsInZ/2+0.5)) * voxSpacingZ)*2;
        rotm[copyNo]  = new G4RotationMatrix(G4ThreeVector(0,1,0), 90*deg);
      }
 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NuLatLightGuideParameterisation::~NuLatLightGuideParameterisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NuLatLightGuideParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const
{
  physVol->SetTranslation(G4ThreeVector(xLightGuide[copyNo],yLightGuide[copyNo],zLightGuide[copyNo]));
  physVol->SetRotation(rotm[copyNo]);
  physVol->SetCopyNo(copyNo);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
