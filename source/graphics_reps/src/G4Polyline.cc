// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Polyline.cc,v 1.4 1999/12/15 14:50:36 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// 
// John Allison  July 1995

#include "G4Polyline.hh"

G4Polyline::G4Polyline (const G4VVisPrim& prim):
  G4VVisPrim (prim)
{}

G4Visible & G4Polyline::operator = (const G4Visible &right) {
  return G4Visible::operator = (right);
}

G4VVisPrim & G4Polyline::operator = (const G4VVisPrim &right) {
  return G4VVisPrim::operator = (right);
}

G4Polyline & G4Polyline::operator = (const G4Polyline &right) {
  if (&right == this) return *this;
  G4VVisPrim::operator = (right);
  return *this;
}

G4std::ostream& operator << (G4std::ostream& os, const G4Polyline& line) {
  os << "G4Polyline: ";
  os << '\n' << (G4VVisPrim) line;
  os << '\n' << (G4Point3DList) line;
  return os;
}