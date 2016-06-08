// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PlaneCreator.hh,v 1.2 2000/01/21 13:45:28 gcosmo Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// ----------------------------------------------------------------------
// Class G4PlaneCreator
//
// Class description:
//
//

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------
#ifndef G4PLANECREATOR_HH
#define G4PLANECREATOR_HH

#include "G4GeometryCreator.hh"

class G4PlaneCreator: private G4GeometryCreator 
{
  public:

  // Constructor & destructor

    G4PlaneCreator();
    ~G4PlaneCreator();

  // Member functions

    void CreateG4Geometry(STEPentity&);
    void CreateSTEPGeometry(void* G4obj);
    G4String Name() { return "Plane"; }
    static G4PlaneCreator GetInstance() { return csc; }

  // Members

  private:

    static G4PlaneCreator csc;
};

#endif