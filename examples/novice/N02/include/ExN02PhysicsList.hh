// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02PhysicsList.hh,v 1.6 2000/02/29 12:01:52 maire Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This class is a class derived from G4VUserPhysicsList
//      for constructing all particles and processes.
//
//	History
//        first version              10  Jan. 1998 by H.Kurashige
// ------------------------------------------------------------
#ifndef ExN02PhysicsList_h
#define ExN02PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class G4PhotoElectricEffect;
class G4ComptonScattering;
class G4GammaConversion;

class G4MultipleScattering;

class G4eIonisation;
class G4eBremsstrahlung;
class G4eplusAnnihilation;

class G4MuIonisation;
class G4MuBremsstrahlung;
class G4MuPairProduction;

class G4hIonisation;

class ExN02PhysicsList: public G4VUserPhysicsList
{
  public:
    ExN02PhysicsList();
    virtual ~ExN02PhysicsList();

  protected:
    // Construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

    // 
    virtual void SetCuts();

   
  protected:
    // these methods Construct particles 
    virtual void ConstructBosons();
    virtual void ConstructLeptons();
    virtual void ConstructMesons();
    virtual void ConstructBaryons();

  protected:
  // these methods Construct physics processes and register them
    virtual void ConstructGeneral();
    virtual void ConstructEM();
    
  private:

    G4PhotoElectricEffect* thePhotoElectricEffect;
    G4ComptonScattering*   theComptonScattering;
    G4GammaConversion*     theGammaConversion;
    
    G4MultipleScattering*  theeminusMultipleScattering;
    G4eIonisation*         theeminusIonisation;
    G4eBremsstrahlung*     theeminusBremsstrahlung;
    
    G4MultipleScattering*  theeplusMultipleScattering;
    G4eIonisation*         theeplusIonisation;
    G4eBremsstrahlung*     theeplusBremsstrahlung;
    G4eplusAnnihilation*   theeplusAnnihilation;
};

#endif


