// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisExtent.hh,v 1.4 1999/12/15 14:50:35 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// 
// A.Walkden 28/11/95

// Class Description:
// G4VisExtent is a class to prototype extent functions for use in
// instantiating Visualisation windows with an appropriate field of
// view for the object(s) being drawn.
// Class Description - End:

#ifndef G4VISEXTENT_HH
#define G4VISEXTENT_HH

#include "globals.hh"
#include "G4Point3D.hh"

class G4VisExtent
{
public: // With description

  G4VisExtent (G4double xmin = 0., G4double xmax = 0., 
               G4double ymin = 0., G4double ymax = 0., 
               G4double zmin = 0., G4double zmax = 0.);
  G4VisExtent (const G4Point3D& centre, G4double radius);
  ~G4VisExtent ();
  G4bool operator != (const G4VisExtent& e) const;
  G4double  GetXmin         () const;
  G4double  GetXmax         () const;
  G4double  GetYmin         () const;
  G4double  GetYmax         () const;
  G4double  GetZmin         () const;
  G4double  GetZmax         () const;
  G4Point3D GetExtentCentre () const;
  G4Point3D GetExtentCenter () const;
  G4double  GetExtentRadius () const;
  void SetXmin (G4double xmin);
  void SetXmax (G4double xmax);
  void SetYmin (G4double ymin);
  void SetYmax (G4double ymax);
  void SetZmin (G4double zmin);
  void SetZmax (G4double zmax);
  friend G4std::ostream& operator << (G4std::ostream& os, const G4VisExtent& e);

private:
  G4double fXmin, fXmax, fYmin, fYmax, fZmin, fZmax;
};

inline G4double G4VisExtent::GetXmin () const { return fXmin; }
inline G4double G4VisExtent::GetXmax () const { return fXmax; }
inline G4double G4VisExtent::GetYmin () const { return fYmin; }
inline G4double G4VisExtent::GetYmax () const { return fYmax; }
inline G4double G4VisExtent::GetZmin () const { return fZmin; }
inline G4double G4VisExtent::GetZmax () const { return fZmax; }

inline G4Point3D G4VisExtent::GetExtentCenter () const {
  return GetExtentCentre ();
}

inline void G4VisExtent::SetXmin (G4double xmin) {fXmin = xmin;}
inline void G4VisExtent::SetXmax (G4double xmax) {fXmax = xmax;}
inline void G4VisExtent::SetYmin (G4double ymin) {fYmin = ymin;}
inline void G4VisExtent::SetYmax (G4double ymax) {fYmax = ymax;}
inline void G4VisExtent::SetZmin (G4double zmin) {fZmin = zmin;}
inline void G4VisExtent::SetZmax (G4double zmax) {fZmax = zmax;}

#endif