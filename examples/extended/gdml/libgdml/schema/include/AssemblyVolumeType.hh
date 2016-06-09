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
// $Id: AssemblyVolumeType.hh,v 1.2 2002/06/03 12:09:33 radoone Exp $
// GEANT4 tag $Name: geant4-05-02-patch-01 $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef ASSEMBLYVOLUMETYPE_H
#define ASSEMBLYVOLUMETYPE_H 1

#include "IdentifiableVolumeType.hh"
#include "SinglePlacementType.hh"
#include "ContentGroup.hh"

class AssemblyVolumeType : public IdentifiableVolumeType {
public:
  AssemblyVolumeType() {
  }
  ~AssemblyVolumeType() {
  }
  
  const ContentSequence* get_content() const {
    return &m_sequence;
  }

  void add_content( const std::string& tag, SAXObject* so ) {
    ContentGroup::ContentItem ci = { tag, so };
    m_sequence.add_content( ci );
  }
private:
  ContentSequence m_sequence;
};

#endif // ASSEMBLYVOLUMETYPE_H