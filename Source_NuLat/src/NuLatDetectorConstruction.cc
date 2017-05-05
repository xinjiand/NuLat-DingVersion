#include "NuLatDetectorConstruction.hh"
#include "NuLatVoxelParameterisation.hh"
#include "NuLatLightGuideParameterisation.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"

#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4IntersectionSolid.hh"
#include "NuLatVoxelSensitiveDetector.hh"
#include "NuLatPhotoCathodeSensitiveDetector.hh"
#include "G4PVParameterised.hh"

#include "G4VSensitiveDetector.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"



using namespace std;


G4VPhysicalVolume* NuLatDetectorConstruction::Construct()
{
  NuLatMaterials = new Materials();

  buildExperimentalHall();

  return experimentalHallPhys;
}


void NuLatDetectorConstruction::buildExperimentalHall()
{
  G4bool checkOverlaps = true;
  G4int copyNumber;

  experimentalHallPhys = new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, 0.0 ),
              WorldVolume(), "AirWorld", 0, false, copyNumber=0, checkOverlaps );
}

/*Create a world logical volume*/
G4LogicalVolume* NuLatDetectorConstruction::WorldVolume()
{
  G4bool checkOverlaps = true;
  photoCathodeIDNum=0;



  G4double totalNumberOfVoxels = nOfVoxelsInX*nOfVoxelsInY*nOfVoxelsInZ;


  G4double NuLatVoxelatedCalorimeterBoxXDimension =
    nOfVoxelsInX * (voxelXDimension+voxelSpacingXDimension)+voxelSpacingXDimension;

  G4double NuLatVoxelatedCalorimeterBoxYDimension =
    nOfVoxelsInY * (voxelYDimension+voxelSpacingYDimension)+voxelSpacingYDimension;

  G4double NuLatVoxelatedCalorimeterBoxZDimension =
    nOfVoxelsInZ * (voxelZDimension+voxelSpacingZDimension)+voxelSpacingZDimension;



  G4Box* world = new G4Box( "AirWorld",
                  NuLatVoxelatedCalorimeterBoxXDimension,
                  NuLatVoxelatedCalorimeterBoxYDimension,
                  NuLatVoxelatedCalorimeterBoxZDimension);
  experimentalHallLog = new G4LogicalVolume( world, NuLatMaterials->air, "AirWorld", 0, 0, 0);


    // NuLat Voxelated Calorimeter
    G4VSolid* NuLatVoxelatedCalorimeterSolid
      = new G4Box("NuLatVoxelatedCalorimeterBox",
                  NuLatVoxelatedCalorimeterBoxXDimension/2,
                  NuLatVoxelatedCalorimeterBoxYDimension/2,
                  NuLatVoxelatedCalorimeterBoxZDimension/2);

    G4LogicalVolume* NuLatVoxelatedCalorimeterLogical
      = new G4LogicalVolume(NuLatVoxelatedCalorimeterSolid,
                            NuLatMaterials->air,
                            "NuLatVoxelatedCalorimeterLogical");
                            
    new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),NuLatVoxelatedCalorimeterLogical,
                      "NuLatVoxelatedCalorimeterPhysical",experimentalHallLog,
                      false,0,checkOverlaps);
    
    // NuLat voxels
    G4VSolid* VoxelSolid
      = new G4Box("Voxel",voxelXDimension/2, voxelYDimension/2, voxelZDimension/2);
      
    VoxelLogical
      = new G4LogicalVolume(VoxelSolid,NuLatMaterials->Li6PVT_5TenthsOfAPercent,"VoxelLogical");//////////
    G4VPVParameterisation* voxelParam =
      new NuLatVoxelParameterisation(nOfVoxelsInX, nOfVoxelsInY, nOfVoxelsInZ,
        voxelXDimension/2, voxelYDimension/2, voxelZDimension/2,
        voxelSpacingXDimension/2, voxelSpacingYDimension/2, voxelSpacingZDimension/2);
      
    new G4PVParameterised("voxelPhysical",VoxelLogical,NuLatVoxelatedCalorimeterLogical,
                          kXAxis,totalNumberOfVoxels,voxelParam);

    // NuLat Light Guides
    LightGuideAndPMTLog = LightGuideAndPMT((voxelXDimension+voxelSpacingXDimension-0.127*cm), 4.6*cm,
                                           (voxelYDimension+voxelSpacingYDimension-0.127*cm), 4.6*cm,
                                            LightGuideTaperLength);
    
    MirrorAndPMTLog = MirrorAndPMT ((voxelXDimension+voxelSpacingXDimension-0.127*cm),
    								   (voxelYDimension+voxelSpacingYDimension-0.127*cm),
									   1.811*2.54/2*cm ,
									   0.1*mm) ;
    G4double MirrorAndPMTLength = 0.1*mm+20*cm ;

    // NuLat Mirror test


    // NuLat LightGuide z+ Bank
    G4VSolid* NuLatLightGuideZBankPlusSolid
      = new G4Box("NuLatLightGuideZBankPlusBox",
                  NuLatVoxelatedCalorimeterBoxXDimension/2,
                  NuLatVoxelatedCalorimeterBoxYDimension/2,
                  (MirrorAndPMTLength)/2);
      
    NuLatLightGuideZBankPlusLogical
      = new G4LogicalVolume(NuLatLightGuideZBankPlusSolid,
                            NuLatMaterials->air,
                            "NuLatLightGuideZBankPlusLogical");
                            
    new G4PVPlacement(0,G4ThreeVector(0.,0.,NuLatVoxelatedCalorimeterBoxZDimension/2+(MirrorAndPMTLength)/2 +.010*2.54*cm),NuLatLightGuideZBankPlusLogical,
                      "NuLatLightGuideZBankPlusPhysical",experimentalHallLog,
                      false,0,checkOverlaps);


    G4VPVParameterisation* lightGuideParamZBankPlus =
      new NuLatLightGuideParameterisation(nOfVoxelsInX, nOfVoxelsInY, nOfVoxelsInZ,
        voxelXDimension/2, voxelYDimension/2, voxelZDimension/2,
        voxelSpacingXDimension/2, voxelSpacingYDimension/2, voxelSpacingZDimension/2, 0, 0, 1);
      
    new G4PVParameterised("LightGuidePhysicalZPlus", MirrorAndPMTLog,NuLatLightGuideZBankPlusLogical,
                          kZAxis,nOfVoxelsInX*nOfVoxelsInY,lightGuideParamZBankPlus);


    // NuLat LightGuide z- Bank
    G4VSolid* NuLatLightGuideZBankMinusSolid
      = new G4Box("NuLatLightGuideZBankMinusBox",
                  NuLatVoxelatedCalorimeterBoxXDimension/2,
                  NuLatVoxelatedCalorimeterBoxYDimension/2,
                  (lightGuideWithPMTLength)/2);
      
    NuLatLightGuideZBankMinusLogical
      = new G4LogicalVolume(NuLatLightGuideZBankMinusSolid,
                            NuLatMaterials->air,
                            "NuLatLightGuideZBankMinusLogical");
                            
    new G4PVPlacement(0,G4ThreeVector(0.,0.,-1*(NuLatVoxelatedCalorimeterBoxZDimension/2+(lightGuideWithPMTLength)/2 +.010*2.54*cm)),NuLatLightGuideZBankMinusLogical,
                      "NuLatLightGuideZBankMinusPhysical",experimentalHallLog,
                      false,0,checkOverlaps);


  G4VPVParameterisation* lightGuideParamZBankMinus =
      new NuLatLightGuideParameterisation(nOfVoxelsInX, nOfVoxelsInY, nOfVoxelsInZ,
        voxelXDimension/2, voxelYDimension/2, voxelZDimension/2,
        voxelSpacingXDimension/2, voxelSpacingYDimension/2, voxelSpacingZDimension/2, 0, 0, -1);
      
    new G4PVParameterised("LightGuidePhysicalZMinus",LightGuideAndPMTLog,NuLatLightGuideZBankMinusLogical,
                          kZAxis,nOfVoxelsInX*nOfVoxelsInY,lightGuideParamZBankMinus);
                          
                          
                          
                          
                          
                          

    // NuLat LightGuide y+ Bank
    G4VSolid* NuLatLightGuideYBankPlusSolid 
      = new G4Box("NuLatLightGuideYBankPlusBox",
                  NuLatVoxelatedCalorimeterBoxXDimension/2,
                  (lightGuideWithPMTLength)/2,
                  NuLatVoxelatedCalorimeterBoxZDimension/2);
      
    NuLatLightGuideYBankPlusLogical
      = new G4LogicalVolume(NuLatLightGuideYBankPlusSolid,
                            NuLatMaterials->air,
                            "NuLatLightGuideYBankPlusLogical");
                            
    new G4PVPlacement(0,G4ThreeVector(0.,NuLatVoxelatedCalorimeterBoxYDimension/2+(lightGuideWithPMTLength)/2 +.010*2.54*cm,0.),NuLatLightGuideYBankPlusLogical,
                      "NuLatLightGuideYBankPlusPhysical",experimentalHallLog,
                      false,0,checkOverlaps);


  G4VPVParameterisation* lightGuideParamYBankPlus =
      new NuLatLightGuideParameterisation(nOfVoxelsInX, nOfVoxelsInY, nOfVoxelsInZ,
        voxelXDimension/2, voxelYDimension/2, voxelZDimension/2,
        voxelSpacingXDimension/2, voxelSpacingYDimension/2, voxelSpacingZDimension/2, 0, 1, 0);
      
    new G4PVParameterised("LightGuidePhysicalZPlus",LightGuideAndPMTLog,NuLatLightGuideYBankPlusLogical,
                          kZAxis,nOfVoxelsInX*nOfVoxelsInZ,lightGuideParamYBankPlus);


    // NuLat LightGuide y- Bank
    G4VSolid* NuLatLightGuideYBankMinusSolid 
      = new G4Box("NuLatLightGuideYBankPlusBox",
                  NuLatVoxelatedCalorimeterBoxXDimension/2,
                  (lightGuideWithPMTLength)/2,
                  NuLatVoxelatedCalorimeterBoxZDimension/2);
      
    NuLatLightGuideYBankMinusLogical
      = new G4LogicalVolume(NuLatLightGuideYBankMinusSolid,
                            NuLatMaterials->air,
                            "NuLatLightGuideYBankMinusLogical");

    new G4PVPlacement(0,G4ThreeVector(0.,-1*(NuLatVoxelatedCalorimeterBoxYDimension/2+(lightGuideWithPMTLength)/2 +.010*2.54*cm),0.),NuLatLightGuideYBankMinusLogical,
                      "NuLatLightGuideYBankMinusPhysical",experimentalHallLog,
                      false,0,checkOverlaps);


  G4VPVParameterisation* lightGuideParamYBankMinus =
      new NuLatLightGuideParameterisation(nOfVoxelsInX, nOfVoxelsInY, nOfVoxelsInZ,
        voxelXDimension/2, voxelYDimension/2, voxelZDimension/2,
        voxelSpacingXDimension/2, voxelSpacingYDimension/2, voxelSpacingZDimension/2, 0, -1, 0);
      
    new G4PVParameterised("LightGuidePhysicalZMinus",LightGuideAndPMTLog,NuLatLightGuideYBankMinusLogical,
                          kZAxis,nOfVoxelsInX*nOfVoxelsInZ,lightGuideParamYBankMinus);
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
    // NuLat LightGuide x+ Bank
    G4VSolid* NuLatLightGuideXBankPlusSolid 
      = new G4Box("NuLatLightGuideXBankPlusBox",
                  (lightGuideWithPMTLength)/2,
                  NuLatVoxelatedCalorimeterBoxYDimension/2,
                  NuLatVoxelatedCalorimeterBoxZDimension/2);
      
    NuLatLightGuideXBankPlusLogical
      = new G4LogicalVolume(NuLatLightGuideXBankPlusSolid,
                            NuLatMaterials->air,
                            "NuLatLightGuideXBankPlusLogical");
                            
    new G4PVPlacement(0,G4ThreeVector(NuLatVoxelatedCalorimeterBoxXDimension/2+(lightGuideWithPMTLength)/2 +.010*2.54*cm,0.,0.),NuLatLightGuideXBankPlusLogical,
                      "NuLatLightGuideXBankPlusPhysical",experimentalHallLog,
                      false,0,checkOverlaps);


  G4VPVParameterisation* lightGuideParamXBankPlus =
      new NuLatLightGuideParameterisation(nOfVoxelsInX, nOfVoxelsInY, nOfVoxelsInZ,
        voxelXDimension/2, voxelYDimension/2, voxelZDimension/2,
        voxelSpacingXDimension/2, voxelSpacingYDimension/2, voxelSpacingZDimension/2, 1, 0, 0);
      
    new G4PVParameterised("LightGuidePhysicalZPlus",LightGuideAndPMTLog,NuLatLightGuideXBankPlusLogical,
                          kZAxis,nOfVoxelsInZ*nOfVoxelsInY,lightGuideParamXBankPlus);


    // NuLat LightGuide x- Bank
    G4VSolid* NuLatLightGuideXBankMinusSolid 
      = new G4Box("NuLatLightGuideXBankPlusBox",
                  (lightGuideWithPMTLength)/2,
                  NuLatVoxelatedCalorimeterBoxYDimension/2,
                  NuLatVoxelatedCalorimeterBoxZDimension/2);
      
    NuLatLightGuideXBankMinusLogical
      = new G4LogicalVolume(NuLatLightGuideXBankMinusSolid,
                            NuLatMaterials->air,
                            "NuLatLightGuideXBankMinusLogical");
                            
    new G4PVPlacement(0,G4ThreeVector(-1*(NuLatVoxelatedCalorimeterBoxXDimension/2+(lightGuideWithPMTLength)/2 +.010*2.54*cm),0.,0.),NuLatLightGuideXBankMinusLogical,
                      "NuLatLightGuideXBankMinusPhysical",experimentalHallLog,
                      false,0,checkOverlaps);


  G4VPVParameterisation* lightGuideParamXBankMinus =
      new NuLatLightGuideParameterisation(nOfVoxelsInX, nOfVoxelsInY, nOfVoxelsInZ,
        voxelXDimension/2, voxelYDimension/2, voxelZDimension/2,
        voxelSpacingXDimension/2, voxelSpacingYDimension/2, voxelSpacingZDimension/2, -1, 0, 0);
      
    new G4PVParameterised("LightGuidePhysicalZMinus",LightGuideAndPMTLog,NuLatLightGuideXBankMinusLogical,
                          kZAxis,nOfVoxelsInZ*nOfVoxelsInY,lightGuideParamXBankMinus);


  // visualization attributes
    G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(0.9,0.9,0.0));
    VoxelLogical->SetVisAttributes(visAttributes);
    NuLatVisAttributes.push_back(visAttributes);
  
  //visAttributes = new G4VisAttributes(false);
    visAttributes = new G4VisAttributes(G4Colour(0.9,0.0,0.0));
    NuLatVoxelatedCalorimeterLogical->SetVisAttributes(visAttributes);
    NuLatVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(false);
    experimentalHallLog->SetVisAttributes(visAttributes);
    NuLatVisAttributes.push_back(visAttributes);


  return(experimentalHallLog);

}



G4LogicalVolume* NuLatDetectorConstruction::HamamatsuR10533()
{


  G4bool checkOverlaps = true;
	//	The PMTs are semi spherical quartz with vacuum inside, 
	//	with a cylindrical mu-metal shield
	G4double pmtRadius=23*mm;
	G4double pmtGlassThickness = 4*mm;
  G4double pmtGlassMinimumThickness = 0.8*mm;
  G4double pmtPhotoCathodeRadius=529.25*mm;
	G4double pmtPhotoCathodeThickness = 1*mm;

	G4double muMetalRadius = 3.05*cm;
	G4double muMetalHeight = 19.863*cm;
	G4double muMetalThickness = 0.05*cm;

// muMetal
	G4Tubs* muMetal_solid = new G4Tubs("muMetal_solid",
	                                    muMetalRadius-muMetalThickness,
	                                    muMetalRadius, 
	                                    muMetalHeight/2., 
	                                    0., 
	                                    360*deg);
  muMetal_Log = new G4LogicalVolume(muMetal_solid,
                                              NuLatMaterials->muMetal,
                                              "muMetal_Log");

// PMT lens
	G4Tubs* pmtGlass_solid = new G4Tubs("pmtGlass_solid", 
	                                     0,
	                                     pmtRadius, 
	                                     pmtGlassThickness/2., 
	                                     0., 
	                                     360*deg);
	                                     
	G4Sphere* pmtConvexSurface_solid = new G4Sphere("pmtConvexSurface_solid",
	                                                0.,
	                                                pmtPhotoCathodeRadius, 
                                                   0.*deg,
	                                                360.*deg,
	                                                0.*deg,
	                                                360.*deg);
	                                                
  G4SubtractionSolid* pmtLens_solid = new G4SubtractionSolid("pmtLens_solid", 
                                                  pmtGlass_solid,
	                                                pmtConvexSurface_solid,
                                                  0,
	                                                G4ThreeVector(0, 0, pmtPhotoCathodeRadius-pmtGlassThickness/2+pmtGlassMinimumThickness));
	                                                
  pmtLens_Log = new G4LogicalVolume(pmtLens_solid,
                                              NuLatMaterials->borosilicateGlass,
                                              "pmtLens_Log");

// photocathode
	G4Sphere* pmtPhotoCathode_solid = new G4Sphere("pmtPhotoCathode_solid",
	                                               pmtPhotoCathodeRadius-pmtPhotoCathodeThickness,
                                                 pmtPhotoCathodeRadius,
	                                               0.*deg,
          	                                     360.*deg,
          	                                     177.51*deg,//177.509
          	                                     180*deg);

  pmtPhotoCathode_Log
                        = new G4LogicalVolume(pmtPhotoCathode_solid,
                                              NuLatMaterials->BeCuPhotoCathode,
                                              "pmtPhotoCathode_Log");

//pmt assembly
  G4Tubs* pmt_solid = new G4Tubs("pmt_solid",
	                                    0,
	                                    muMetalRadius, 
	                                    muMetalHeight/2., 
	                                    0., 
	                                    360*deg);

  G4LogicalVolume * pmt_Log
      = new G4LogicalVolume(pmt_solid,
                            NuLatMaterials->vacuum,
                            "pmt_log");

  new G4PVPlacement(0,G4ThreeVector(0.,0.,0),muMetal_Log,
                      "muMetal",pmt_Log,
                      false,0,checkOverlaps);

  new G4PVPlacement(0,G4ThreeVector(0.,0.,-1*(muMetalHeight/2.- pmtGlassThickness/2)),pmtLens_Log,
                      "pmtLens",pmt_Log,
                      false,0,checkOverlaps);

  new G4PVPlacement(0,
                    G4ThreeVector(0.,0.,pmtPhotoCathodeRadius-muMetalHeight/2+pmtGlassMinimumThickness),//-1*(.- pmtGlassThickness/2)-)
                    pmtPhotoCathode_Log,
                      "pmtPhotoCathode",pmt_Log,
                      false,0,checkOverlaps);

  G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(0,255,0,.3));
  visAttributes->SetForceWireframe(true);
//  visAttributes->SetForceSolid(true);
  muMetal_Log->SetVisAttributes(visAttributes);
  NuLatVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(255,0,0,.3));
//  visAttributes->SetForceWireframe(true);
  visAttributes->SetForceSolid(true);
  pmtLens_Log->SetVisAttributes(visAttributes);
  NuLatVisAttributes.push_back(visAttributes);
  
//  visAttributes = new G4VisAttributes(G4Colour(0,0,255,.3));
//  visAttributes->SetForceWireframe(true);
//  visAttributes->SetForceSolid(true);
//  pmtFace_Log->SetVisAttributes(visAttributes);
//  NuLatVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0,255,0,.3));
  visAttributes->SetForceSolid(true);
  pmtPhotoCathode_Log->SetVisAttributes(visAttributes);
  NuLatVisAttributes.push_back(visAttributes);

  
  
  
  return (pmt_Log);
}


G4LogicalVolume* NuLatDetectorConstruction::LightGuideAndPMT(G4double dx1, G4double dx2, G4double dy1, G4double dy2, G4double dz)
{

  G4VSolid* lightGuideTrd = new G4Trd("LightGuideTrd", dx1/2, dx2/2, dy1/2, dy2/2, dz/2);
  
  G4double r1max = 3.465*2.54/2*cm;
  G4double r1min = 0.0*cm;
  G4double r2max = 1.811*2.54/2*cm;
  G4double r2min = 0.0*cm;
  G4VSolid* lightGuideCone = new G4Cons("LightGuideCone", r1min, r1max, r2min, r2max, dz/2, 0, 360*deg);
  G4bool checkOverlaps = true;
  
  G4Box* lightGuideSquare = new G4Box( "LightGuideSquare", dx1/2, dy1/2, 0.5*cm/2);

  G4IntersectionSolid* lightGuideTrdIntersLightGuideCone = new G4IntersectionSolid("Trd Intersect Cone", lightGuideTrd, lightGuideCone);

    G4VSolid* lightGuideBox 
      = new G4Box("lightGuideBox",dx1/2, dy1/2, (dz+0.5*cm+19.863*cm)/2);

    G4LogicalVolume* lightGuide
      = new G4LogicalVolume(lightGuideBox,
                            NuLatMaterials->air,
                            "lightGuide");
    lightGuideTrdIntersLightGuideConeLog
      = new G4LogicalVolume(lightGuideTrdIntersLightGuideCone,
                            NuLatMaterials->acrylic,
                            "lightGuideTrdIntersLightGuideConeLog");
    lightGuideSquareLog
      = new G4LogicalVolume(lightGuideSquare,
                            NuLatMaterials->acrylic,
                            "lightGuideSquareLog");

    new G4PVPlacement(0,G4ThreeVector(0.,0.,0.5*cm/2 - 19.863*cm/2),lightGuideTrdIntersLightGuideConeLog,
                      "guide001",lightGuide,
                      false,0,checkOverlaps);
   new G4PVPlacement(0,G4ThreeVector(0.,0.,-1.0*dz/2 - 19.863*cm/2),lightGuideSquareLog,
                      "guide002",lightGuide,
                      false,0,checkOverlaps);

   PMTLog = HamamatsuR10533();
   new G4PVPlacement(0,G4ThreeVector(0.,0.,(dz+0.5*cm)/2),PMTLog,
                      "PMT",lightGuide,
                      false,0,checkOverlaps);

  G4VisAttributes* visAttributes = new G4VisAttributes(false);
    lightGuide->SetVisAttributes(visAttributes);
    NuLatVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour(255,0,255,.3));
    lightGuideTrdIntersLightGuideConeLog->SetVisAttributes(visAttributes);

    lightGuideSquareLog->SetVisAttributes(visAttributes);
    NuLatVisAttributes.push_back(visAttributes);

    return(lightGuide);
}

G4LogicalVolume* NuLatDetectorConstruction::MirrorAndPMT(G4double xd, G4double yd, G4double rd, G4double zd)
{
	  //G4VSolid* lightGuideTrd = new G4Trd("LightGuideTrd", dx1/2, dx2/2, dy1/2, dy2/2, dz/2);

	 /* G4double r1max = 3.465*2.54/2*cm;
	  G4double r1min = 0.0*cm;
	  G4double r2max = 1.811*2.54/2*cm;
	  G4double r2min = 0.0*cm;
	  G4VSolid* lightGuideCone = new G4Cons("LightGuideCone", r1min, r1max, r2min, r2max, dz/2, 0, 360*deg);
	  G4bool checkOverlaps = true;

	  G4Box* lightGuideSquare = new G4Box( "LightGuideSquare", dx1/2, dy1/2, 0.5*cm/2); */

	    G4Cons* MirrorCircle = new G4Cons ("MirrorCircle",0*cm, rd , 0*cm, rd, zd/2 , 0, 360*deg);
	    G4Box* MirrorSquare = new G4Box ("MirrorBox" , xd/2 , yd/2 , zd/2);
	    G4SubtractionSolid* MirrorOpen= new G4SubtractionSolid ("MirrorSquare-Circle" , MirrorSquare, MirrorCircle);
	    G4bool checkOverlaps = true ;
	  //G4IntersectionSolid* lightGuideTrdIntersLightGuideCone = new G4IntersectionSolid("Trd Intersect Cone", lightGuideTrd, lightGuideCone);

	    G4VSolid* lightGuideBox
	      = new G4Box("lightGuideBox",xd/2, yd/2, (zd+19.863*cm)/2);

	    G4LogicalVolume* lightGuide
	      = new G4LogicalVolume(lightGuideBox,
	                            NuLatMaterials->air,
	                            "lightGuide");
	    MirrorLogical
	      = new G4LogicalVolume(MirrorSquare,
	                            NuLatMaterials->muMetal,
	                            "MirrorLogical");
	   /* lightGuideSquareLog
	     = new G4LogicalVolume(lightGuideSquare,
	                            NuLatMaterials->acrylic,
	                            "lightGuideSquareLog"); */

	    new G4PVPlacement(0,G4ThreeVector(0.,0.,- 19.863*cm/2),MirrorLogical,
	                      "guide001",lightGuide,
	                      false,0,checkOverlaps);
	   // new G4PVPlacement(0,G4ThreeVector(0.,0.,1.0*dz/2 - 19.863*cm/2),lightGuideSquareLog,
	   //                   "guide002",lightGuide,
	   //                   false,0,checkOverlaps);

	   PMTLog = HamamatsuR10533();
	   new G4PVPlacement(0,G4ThreeVector(0.,0.,zd/2),PMTLog,
	                      "PMT",lightGuide,
	                      false,0,checkOverlaps);

	  G4VisAttributes* visAttributes = new G4VisAttributes(false);
	    lightGuide->SetVisAttributes(visAttributes);
	    NuLatVisAttributes.push_back(visAttributes);

	  visAttributes = new G4VisAttributes(G4Colour(255,0,255,.3));
	    lightGuideTrdIntersLightGuideConeLog->SetVisAttributes(visAttributes);

	    lightGuideSquareLog->SetVisAttributes(visAttributes);
	    NuLatVisAttributes.push_back(visAttributes);

	    return(lightGuide);
}


void NuLatDetectorConstruction::ConstructSDandField()
{
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4String SDname;
    
    G4VSensitiveDetector* nuLatVoxel
      = new NuLatVoxelSensitiveDetector(SDname="/NuLatVoxel",
        nOfVoxelsInX, nOfVoxelsInY, nOfVoxelsInZ,
        voxelXDimension/2, voxelYDimension/2, voxelZDimension/2,
        voxelSpacingXDimension/2, voxelSpacingYDimension/2, voxelSpacingZDimension/2);
    SDman->AddNewDetector(nuLatVoxel);
    VoxelLogical->SetSensitiveDetector(nuLatVoxel);

    G4VSensitiveDetector* nuLatPhotoCathode
      = new NuLatPhotoCathodeSensitiveDetector(SDname="/nuLatPhotoCathode",
        nOfVoxelsInX, nOfVoxelsInY, nOfVoxelsInZ,
        voxelXDimension/2, voxelYDimension/2, voxelZDimension/2,
        voxelSpacingXDimension/2, voxelSpacingYDimension/2, voxelSpacingZDimension/2);
    SDman->AddNewDetector(nuLatPhotoCathode);
    pmtPhotoCathode_Log->SetSensitiveDetector(nuLatPhotoCathode);
}





