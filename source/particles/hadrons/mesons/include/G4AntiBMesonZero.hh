// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AntiBMesonZero.hh,v 1.3 1999/12/15 14:51:02 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//
//      Created,             Hisaya Kurashige, 15 June 1997
// **********************************************************************
//  Change both methods to get the pointer into non-inlined H.Kurashige 4 Aug. 1998
// ----------------------------------------------------------------

// Each class inheriting from G4VMeson
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.

#ifndef G4AntiBMesonZero_h
#define G4AntiBMesonZero_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VMeson.hh"

// ######################################################################
// ###                         AntiBMesonZero                         ###
// ######################################################################

class G4AntiBMesonZero : public G4VMeson
{
 private:
   static G4AntiBMesonZero theAntiBMesonZero;
   static G4double  theAntiBMesonZeroLengthCut;
   static G4double* theAntiBMesonZeroKineticEnergyCuts;

 private: // constructors are hide as private  
   G4AntiBMesonZero(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable
   );

 public:
   virtual ~G4AntiBMesonZero(){}
 
   static G4AntiBMesonZero* AntiBMesonZeroDefinition();
   static G4AntiBMesonZero* AntiBMesonZero();
   static G4double GetCuts() {return theAntiBMesonZeroLengthCut;}   
   static G4double* GetCutsInEnergy() {return theAntiBMesonZeroKineticEnergyCuts;};

   virtual void SetCuts(G4double aCut); 
};


#endif