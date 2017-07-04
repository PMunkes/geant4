//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: pyG4DecayTable.cc 66892 2013-01-17 10:57:59Z gunter $
// ====================================================================
//   pyG4DecayTable.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4DecayTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4DecayTable {

void (G4DecayTable::*f2_Insert)(G4PhaseSpaceDecayChannel*)
  = &G4DecayTable::Insert;
void (G4DecayTable::*f1_Insert)(G4VDecayChannel*)
  = &G4DecayTable::Insert;

}
using namespace pyG4DecayTable; 
// ====================================================================
// module definition
// ====================================================================
void export_G4DecayTable()
{
  class_<G4DecayTable, G4DecayTable*, boost::noncopyable>
    ("G4DecayTable", "decay table")
     // ---
     .def("DumpInfo",   &G4DecayTable::DumpInfo)
     .def("Insert",     f1_Insert)
     .def("Insert",     f2_Insert)
     .def("Entries",    &G4DecayTable::entries)
     ;

     // reduced functionality...
     // ...

}

