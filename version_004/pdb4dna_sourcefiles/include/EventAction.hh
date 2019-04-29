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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file EventAction.hh
/// \brief Definition of the EventAction class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <map>

using namespace std;

class EventActionMessenger;

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction();

public:
  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);

  void AddEdepEvent(G4double edep, G4double partx, G4double party, G4double partz)
  {
    if(fTotalEnergyDeposit == 0 && flag == 0)
    {
      flag = 1;

      fTotalEnergyDeposit += edep;
      xfirst = partx;
      yfirst = party;
      zfirst = partz;
      eatfirst = edep;
      interacter.open("interaction_counter.txt", fstream::app);
      interacter << "1" << G4endl;
      interacter.close();
    }
    else
    {
      fTotalEnergyDeposit += edep;
    }
  };
  G4double GetEdepEvent()
  {
    return fTotalEnergyDeposit;
  };

  void AddEdepToNucleotide(G4int numStrand,G4int numNucl,G4double edep)
  {
    if(numStrand==1)
    {
      fEdepStrand1[numNucl]+=edep;
    }
    else{
      fEdepStrand2[numNucl]+=edep;
    }
  }

  void SetEnergyThresForSSB(G4double val)
  {
    fThresEdepForSSB=val;
  };
  void SetDistanceThresForDSB(G4int val)
  {fThresDistForDSB=val;
  };

private:
  // total energy deposit per event
  G4double fTotalEnergyDeposit;
  // map: first strand (G4int : nucleotide ID, G4double : energy deposit)
  std::map<G4int,G4double>  fEdepStrand1;
  // map: second strand (G4int : nucleotide ID, G4double : energy deposit)
  std::map<G4int,G4double>  fEdepStrand2;
  // Min energy to consider single strand break
  G4double fThresEdepForSSB;
  // Max distance to consider double strand break
  G4int fThresDistForDSB;
  //Flag for the identification of an already going process, for lack of better understanding - Daniel
  G4int flag;
  //Coordinates for the first particle interaction on the event - Daniel
  G4double xfirst;
  G4double yfirst;
  G4double zfirst;
  //Energy deposited in the first electron interaction
  G4double eatfirst;
  //Coordinates for the atom comparison - Daniel
  G4double atomfx;
  G4double atomfy;
  G4double atomfz;
  //Smallest Distances - Daniel
  G4double Odistance;
  G4double Ocomparison;
  G4double NOdistance;
  G4double NOcomparison;
  //Strings for the identification of the atom interacting with the electron - Daniel
  string atomline = "";
  G4double atomtype = 0; //0 is no interaction, 1 is oxigen, 2 is other atom
  G4double closest_o = 0, closest_no = 0, mol_inter = 0; // definition of closest molecule type to oxigen and not oxigen interaction 1 is collagen 2 is elastin
  //Files pointers
  std::ofstream file;
  std::ofstream interacter;
  std::ifstream pdeebfile;



  EventActionMessenger*     fpEventMessenger;

  // Compute Strand breaks from energy deposits in DNA strands
  void ComputeStrandBreaks(G4int*);

  //! return distance between two 3D points
  G4double EADistanceTwo3Dpoints(G4double xA,G4double xB,
      G4double yA,G4double yB,
      G4double zA,G4double zB);


};

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
