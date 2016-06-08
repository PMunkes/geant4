// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateXmViewer.hh,v 1.3 1999/12/15 14:54:04 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// 
// Andrew Walkden  10th February 1997
// G4OpenGLXmViewer : Class derived from G4OpenGLXViewer, to provide
//                  (Motif) widget OpenGL functionality for GEANT4.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLIMMEDIATEXMVIEWER_HH
#define G4OPENGLIMMEDIATEXMVIEWER_HH

#include "G4VViewer.hh"
#include "G4OpenGLImmediateViewer.hh"
#include "G4OpenGLXmViewer.hh"

#include "globals.hh"
#include "g4rw/tvordvec.h"

class G4OpenGLImmediateSceneHandler;

class G4OpenGLImmediateXmViewer:
public G4OpenGLXmViewer, public G4OpenGLImmediateViewer{
  
public:
  G4OpenGLImmediateXmViewer (G4OpenGLImmediateSceneHandler& scene,
			   const G4String& name = "");
  void DrawView ();

};

#endif

#endif