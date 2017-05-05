#include "Randomize.hh"
#include "G4ios.hh"
#include "math.h"

#include "NuLatPrimaryParticleGenerator.hh"


#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryVertex.hh"
//#include "G4GenericMessenger.hh" to be implemented
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>
using namespace std;

NuLatPrimaryParticleGenerator::NuLatPrimaryParticleGenerator() 
{
  G4int n_particle = 1;
  NuLatGenerator = new NuLatPrimaryGeneratorMessenger(this);
  fParticleGun  = new G4ParticleGun(n_particle);
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  fGeantino = particleTable->FindParticle(particleName="geantino");
  fPositron = particleTable->FindParticle(particleName="e+");
  fMuon = particleTable->FindParticle(particleName="mu+");
  fPion = particleTable->FindParticle(particleName="pi+");
  fKaon = particleTable->FindParticle(particleName="kaon+");
  fProton = particleTable->FindParticle(particleName="proton");
  fElectron = particleTable->FindParticle(particleName="e-");
  fGamma = particleTable->FindParticle(particleName="gamma");
  fNeutron = particleTable->FindParticle(particleName="neutron");
}

NuLatPrimaryParticleGenerator::~NuLatPrimaryParticleGenerator()
{
  delete fParticleGun;
//  delete fMessenger; to be implemented
}

void NuLatPrimaryParticleGenerator::GeneratePrimaries( G4Event* event )
{

  if ( EventTypeFlag == "electronTest")
  {
	  GenerateTestParticle( event, fElectron );
  }
  else if ( EventTypeFlag == "positronTest")
  {
	  GenerateTestParticle( event, fPositron );
  }
  else  if ( EventTypeFlag == "IBDTest")
  {
	  GenerateIBDEvent( event );
  }
  else  //default
  {

	  GenerateTestParticle( event, fGeantino );

  }

  fParticleGun->GeneratePrimaryVertex(event);

}






void NuLatPrimaryParticleGenerator::GenerateTestParticle( G4Event* event, G4ParticleDefinition* particle )
{

	G4double Ekin=500*keV;
	G4double posX = 0*cm,  posY = 0*cm,    posZ = 0*cm; //position
    G4double dirX = 0, dirY = 0, dirZ = 1; //Momentum unit vector
    G4double particleTime=0.0*ns;
    G4PrimaryVertex* fvertex= new G4PrimaryVertex();
    G4PrimaryParticle* pp = new G4PrimaryParticle();

    pp->SetParticleDefinition(particle);
    pp->SetMomentumDirection(G4ThreeVector(dirX, dirY, dirZ));
    pp->SetKineticEnergy(Ekin);

    fvertex->SetPrimary(pp);
    fvertex->SetPosition(posX,posY,posZ);
    fvertex->SetT0(0.0);
    fvertex->SetT0(particleTime);

    event->AddPrimaryVertex(fvertex);
}

void NuLatPrimaryParticleGenerator::GenerateIBDEvent( G4Event* event )
{
	G4double Ekin;
	G4double posX = 0*cm,  posY = 0*cm,    posZ = 0*cm; //position
    G4double dirX = 0, dirY = 0, dirZ = 1; //Momentum unit vector
    G4double particleTime=0.0*ns;
    G4PrimaryVertex* fvertex1= new G4PrimaryVertex();
    G4PrimaryVertex* fvertex2= new G4PrimaryVertex();
    G4PrimaryParticle* pp = new G4PrimaryParticle();
    G4PrimaryParticle* pp2 = new G4PrimaryParticle();

    Ekin=GetRndReactorPositronEnergy();

    pp->SetParticleDefinition(fPositron);
    pp->SetMomentumDirection(G4ThreeVector(dirX, dirY, dirZ));
    pp->SetKineticEnergy(Ekin);

    fvertex1->SetPosition(posX,posY,posZ);
    fvertex1->SetT0(0.0);
    fvertex1->SetPrimary(pp);
    fvertex1->SetT0(particleTime);

    event->AddPrimaryVertex(fvertex1);

    Ekin=5*eV;
    pp2->SetParticleDefinition(fNeutron);
    pp2->SetMomentumDirection(G4ThreeVector(dirX, dirY, dirZ));
    pp2->SetKineticEnergy(Ekin);

    fvertex2->SetPosition(posX,posY,posZ);
    fvertex2->SetT0(0.0);
    fvertex2->SetPrimary(pp2);
    fvertex2->SetT0(particleTime);

    event->AddPrimaryVertex(fvertex2);
}

/*******************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
*******************************************************************************/

G4double NuLatPrimaryParticleGenerator::GetRndReactorPositronEnergy()
{
  //Reference journals.aps.org/prd/pdf/10.1103/PhysRevD.34.2621 -- Table IV Y10(Ee+)
  G4double ReactorFlux[26]=
      {  0.000,
         0.150,   0.500,    0.970,    1.175,    1.302,
         1.328,   1.281,    1.179,    1.020,    0.867,
         0.744,   0.617,    0.509,    0.406,    0.307,
         0.239,   0.174,    0.110,    0.050,    0.025,
         0.012,   0.006,    0.003,    0.001,    0.000
      }; //counts/d  first and last data points are extrapolated from reference
  G4double ReactorFluxEnergy[26]=
      { 0.000*MeV,
        0.267*MeV,0.572*MeV,0.877*MeV,1.182*MeV,1.487*MeV,
        1.792*MeV,2.097*MeV,2.402*MeV,2.707*MeV,3.012*MeV,
        3.317*MeV,3.622*MeV,3.927*MeV,4.232*MeV,4.537*MeV,
        4.842*MeV,5.147*MeV,5.452*MeV,5.757*MeV,6.062*MeV,
        6.367*MeV,6.672*MeV,6.977*MeV,7.282*MeV,7.587*MeV
      }; //first and last data points are extrapolated from reference

  G4double EkinOfPositron=-1*MeV;
  G4double rndReactorPositronEnergy,
           leftOfRndReactorPositronEnergy,
           rightOfRndReactorPositronEnergy;
  G4double rndReactorFlux;
  G4double ReactorFluxAtRNDEnergy;
  G4int i;



  while (EkinOfPositron<0.0*MeV)
  {
    rndReactorPositronEnergy = (G4UniformRand() * (ReactorFluxEnergy[25]-ReactorFluxEnergy[0]) + ReactorFluxEnergy[0] );
    rndReactorFlux = (G4UniformRand()*1.328);

    if(rndReactorPositronEnergy==ReactorFluxEnergy[0] || rndReactorPositronEnergy==ReactorFluxEnergy[25])
    {
      EkinOfPositron = -1*MeV; // flux is 0 at end points
    }
    else
    {
      i=0;
      leftOfRndReactorPositronEnergy = ReactorFluxEnergy[i];
      rightOfRndReactorPositronEnergy = ReactorFluxEnergy[i+1];

      while(!(rndReactorPositronEnergy > leftOfRndReactorPositronEnergy &&
        rndReactorPositronEnergy < rightOfRndReactorPositronEnergy))
      {
        i++;
        leftOfRndReactorPositronEnergy = ReactorFluxEnergy[i];
        rightOfRndReactorPositronEnergy = ReactorFluxEnergy[i+1];
      }

      G4double slope =( ReactorFlux[i] - ReactorFlux[i+1] ) / ( ReactorFluxEnergy[i] - ReactorFluxEnergy[i+1] );
      G4double yIntercept = ReactorFlux[i+1] - slope * ReactorFluxEnergy[i+1];
      ReactorFluxAtRNDEnergy = slope*rndReactorPositronEnergy + yIntercept;

      if(rndReactorFlux < ReactorFluxAtRNDEnergy )EkinOfPositron=rndReactorPositronEnergy;
    }
  }
  return(EkinOfPositron);
}



/*******************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
*******************************************************************************/

G4double NuLatPrimaryParticleGenerator::GetRndReactorBackgroundGammaEnergyInMiniTimeCubeCaveAtNIST()
{
  //Reference journals.aps.org/prd/pdf/10.1103/PhysRevD.34.2621
  G4int ReactorFlux[81]=
      {
      100000,  70000,  45000,  32500,  26250,  23125,  21563,  20781,  20391,  20195,
       17509,  17220,  16959,  16714,  16478,  16248,  16023,  15802,  15584,  15369, 
       15157,  14949,  14743,  14540,  14340,  14142,  13947,  13755,  13566,  13379, 
       13195,  13013,  12834,  12658,  12483,  12311,  12142,  11975,  11810,  11647, 
       11487,  11329,  11173,  11019,  10867,  10718,  10570,  10425,  11309,  12168, 
       13000,  13807,  14590,  15348,  16083,  16795,  17484,  18150,  18796,  19419,
       20023, 200227, 151743,  11500,   8715,   6605,   5006,   3794,   2875,   2179, 
        1651,   1251,    948,    719,    545,   4128,    313,    237,    180,    136, 
           0  
      }; // flux/cm^2/s  first and last data points are extrapolated from reference
  G4double ReactorFluxEnergy[81]=
      {         
      0.00*MeV, 0.10*MeV, 0.20*MeV, 0.30*MeV, 0.40*MeV, 0.50*MeV, 0.60*MeV, 0.70*MeV, 0.80*MeV, 0.90*MeV, 
      1.00*MeV, 1.10*MeV, 1.20*MeV, 1.30*MeV, 1.40*MeV, 1.50*MeV, 1.60*MeV, 1.70*MeV, 1.80*MeV, 1.90*MeV, 
      2.00*MeV, 2.10*MeV, 2.20*MeV, 2.30*MeV, 2.40*MeV, 2.50*MeV, 2.60*MeV, 2.70*MeV, 2.80*MeV, 2.90*MeV, 
      3.00*MeV, 3.10*MeV, 3.20*MeV, 3.30*MeV, 3.40*MeV, 3.50*MeV, 3.60*MeV, 3.70*MeV, 3.80*MeV, 3.90*MeV, 
      4.00*MeV, 4.10*MeV, 4.20*MeV, 4.30*MeV, 4.40*MeV, 4.50*MeV, 4.60*MeV, 4.70*MeV, 4.80*MeV, 4.90*MeV, 
      5.00*MeV, 5.10*MeV, 5.20*MeV, 5.30*MeV, 5.40*MeV, 5.50*MeV, 5.60*MeV, 5.70*MeV, 5.80*MeV, 5.90*MeV, 
      6.00*MeV, 6.10*MeV, 6.20*MeV, 6.30*MeV, 6.40*MeV, 6.50*MeV, 6.60*MeV, 6.70*MeV, 6.80*MeV, 6.90*MeV, 
      7.00*MeV, 7.10*MeV, 7.20*MeV, 7.30*MeV, 7.40*MeV, 7.50*MeV, 7.60*MeV, 7.70*MeV, 7.80*MeV, 7.90*MeV, 
      10.00*MeV
      };

  G4double EkinOfGamma=-1*MeV;
  G4double  rndReactorGammaEnergy,
            leftOfRndReactorGammaEnergy,
            rightOfRndReactorGammaEnergy;
  G4double rndReactorGammaFlux;
  G4double reactorGammaFluxAtRNDEnergy;
  G4int i;



  while (EkinOfGamma<0.0*MeV)
  {
    rndReactorGammaEnergy = (G4UniformRand() * (ReactorFluxEnergy[80]-ReactorFluxEnergy[0]) + ReactorFluxEnergy[0] );
    rndReactorGammaFlux = (G4UniformRand()*200227);

    if(rndReactorGammaEnergy==ReactorFluxEnergy[80])
    {
      EkinOfGamma = -1*MeV; // flux is 0 at end points
    }
    else
    {
      i=0;
      leftOfRndReactorGammaEnergy = ReactorFluxEnergy[i];
      rightOfRndReactorGammaEnergy = ReactorFluxEnergy[i+1];

      while(!(rndReactorGammaEnergy > leftOfRndReactorGammaEnergy &&
        rndReactorGammaEnergy < rightOfRndReactorGammaEnergy))
      {
        i++;
        leftOfRndReactorGammaEnergy = ReactorFluxEnergy[i];
        rightOfRndReactorGammaEnergy = ReactorFluxEnergy[i+1];
      }

      G4double slope =( ReactorFlux[i] - ReactorFlux[i+1] ) / ( ReactorFluxEnergy[i] - ReactorFluxEnergy[i+1] );
      G4double yIntercept = ReactorFlux[i+1] - slope * ReactorFluxEnergy[i+1];
      reactorGammaFluxAtRNDEnergy = slope*rndReactorGammaEnergy + yIntercept;

      if(rndReactorGammaFlux < reactorGammaFluxAtRNDEnergy )EkinOfGamma=rndReactorGammaEnergy;
    }
  }
  return(EkinOfGamma);
}



/*******************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
*******************************************************************************/

G4double NuLatPrimaryParticleGenerator::GetRndSpecialNuclearMaterialNeutronFissionEnergy()
{
  //Reference P239 spectrum Eq 13 of Reference http://www.physicsegypt.org/jnrp/v6n003.pdf
  G4int specialNuclearMaterialNeutronFissionFlux[29]=
      {
      0.021295480617, 0.042500438503, 0.056103299644, 0.066913747572, 
      0.131010179356, 0.169661677782, 0.198515503359, 0.320938149138, 
      0.343192103628, 0.331577601354, 0.230722480674, 0.139035361745, 
      0.078992228948, 0.043453934388, 0.023421216922, 0.012447226156, 
      0.006547240332, 0.003416837095, 0.001772118453, 0.000914489703, 
      0.000469962093, 0.000240676458, 0.000122889756, 0.000062587416, 
      0.000031804686, 0.000016130407, 0.000008166704, 0.000004128356, 
      0.000002084035 
      }; 
  G4double specialNuclearMaterialNeutronFissionEnergy[29]=
      {         
      0.001*MeV, 0.004*MeV, 0.007*MeV, 0.010*MeV, 
      0.040*MeV, 0.070*MeV, 0.100*MeV, 0.400*MeV, 
      0.700*MeV, 1.000*MeV, 2.000*MeV, 3.000*MeV, 
      4.000*MeV, 5.000*MeV, 6.000*MeV, 7.000*MeV, 
      8.000*MeV, 9.000*MeV, 10.00*MeV, 11.00*MeV, 
      12.00*MeV, 13.00*MeV, 14.00*MeV, 15.00*MeV, 
      16.00*MeV, 17.00*MeV, 18.00*MeV, 19.00*MeV, 
      20.00*MeV
      };

  G4double EkinOfNeutron=-1*MeV;
  G4double  rndSpecialNuclearMaterialNeutronFissionEnergy,
            leftOfRndSpecialNuclearMaterialNeutronFissionEnergy,
            rightOfRndSpecialNuclearMaterialNeutronFissionEnergy;
  G4double rndSpecialNuclearMaterialNeutronFissionFlux;
  G4double neutronFluxAtRNDEnergy;
  G4int i;



  while (EkinOfNeutron<0.0*MeV)
  {
    rndSpecialNuclearMaterialNeutronFissionEnergy = (G4UniformRand() * (specialNuclearMaterialNeutronFissionEnergy[27]-specialNuclearMaterialNeutronFissionEnergy[0]) + specialNuclearMaterialNeutronFissionEnergy[0] );
    rndSpecialNuclearMaterialNeutronFissionFlux = (G4UniformRand()*0.343192103628);


    i=0;
    leftOfRndSpecialNuclearMaterialNeutronFissionEnergy = specialNuclearMaterialNeutronFissionEnergy[i];
    rightOfRndSpecialNuclearMaterialNeutronFissionEnergy = specialNuclearMaterialNeutronFissionEnergy[i+1];

    while(!(rndSpecialNuclearMaterialNeutronFissionEnergy > leftOfRndSpecialNuclearMaterialNeutronFissionEnergy &&
      rndSpecialNuclearMaterialNeutronFissionEnergy < rightOfRndSpecialNuclearMaterialNeutronFissionEnergy))
    {
      i++;
      leftOfRndSpecialNuclearMaterialNeutronFissionEnergy = specialNuclearMaterialNeutronFissionEnergy[i];
      rightOfRndSpecialNuclearMaterialNeutronFissionEnergy = specialNuclearMaterialNeutronFissionEnergy[i+1];
    }

    G4double slope =( specialNuclearMaterialNeutronFissionFlux[i] - specialNuclearMaterialNeutronFissionFlux[i+1] ) / ( specialNuclearMaterialNeutronFissionEnergy[i] - specialNuclearMaterialNeutronFissionEnergy[i+1] );
    G4double yIntercept = specialNuclearMaterialNeutronFissionFlux[i+1] - slope * specialNuclearMaterialNeutronFissionEnergy[i+1];
    neutronFluxAtRNDEnergy = slope*rndSpecialNuclearMaterialNeutronFissionEnergy + yIntercept;

    if(rndSpecialNuclearMaterialNeutronFissionFlux < neutronFluxAtRNDEnergy )EkinOfNeutron=rndSpecialNuclearMaterialNeutronFissionEnergy;
  }
  return(EkinOfNeutron);
}
