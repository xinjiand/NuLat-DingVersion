#include "NuLatActionInitialization.hh"
#include "NuLatPrimaryParticleGenerator.hh"
#include "NuLatRunAction.hh"
#include "NuLatEventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
NuLatActionInitialization::NuLatActionInitialization(int xVoxels, int yVoxels, int zVoxels)
 : G4VUserActionInitialization()
{
  localNOfVoxelsInX = xVoxels;
  localNOfVoxelsInY = yVoxels;
  localNOfVoxelsInZ = zVoxels;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NuLatActionInitialization::~NuLatActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NuLatActionInitialization::BuildForMaster() const
{
  SetUserAction(new NuLatRunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NuLatActionInitialization::Build() const
{
  NuLatPrimaryParticleGenerator* NuLatPrimaryParticleGenerator001 = new NuLatPrimaryParticleGenerator();
  NuLatRunAction* NuLatRunAction001 = new NuLatRunAction();
  NuLatEventAction* NuLatEventAction001 = new NuLatEventAction(localNOfVoxelsInX, localNOfVoxelsInY, localNOfVoxelsInZ);
  
  SetUserAction(NuLatPrimaryParticleGenerator001);
  SetUserAction(NuLatRunAction001);
  SetUserAction(NuLatEventAction001);
  NuLatEventAction001->NuLatRunActionLocal=NuLatRunAction001;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
