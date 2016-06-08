//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4XXXViewer.cc,v 1.4 2001/11/16 10:50:27 johna Exp $
// GEANT4 tag $Name: geant4-04-00 $

#include "G4XXXViewer.hh"

#include "G4ios.hh"
#include "g4std/strstream"

#include "G4VSceneHandler.hh"

G4XXXViewer::G4XXXViewer
(G4VSceneHandler& sceneHandler, const G4String& name):
  G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name) {}

G4XXXViewer::~G4XXXViewer() {}

void G4XXXViewer::SetView() {
#ifdef G4XXXDEBUG
  G4cout << "G4XXXViewer::SetView() called." << G4endl;
#endif
}

void G4XXXViewer::ClearView() {
#ifdef G4XXXDEBUG
  G4cout << "G4XXXViewer::ClearView() called." << G4endl;
#endif
}

void G4XXXViewer::DrawView() {
#ifdef G4XXXDEBUG
  G4cout << "G4XXXViewer::DrawView() called." << G4endl;
#endif
  NeedKernelVisit ();  // Always need to visit G4 kernel.
  ProcessView ();
}