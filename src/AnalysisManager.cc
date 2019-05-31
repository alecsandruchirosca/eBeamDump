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
// Authors: Susanna Guatelli, susanna@uow.edu.au,
// Authors: Jeremy Davis, jad028@uowmail.edu.au
//

#include <stdlib.h>
#include "AnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

AnalysisManager::AnalysisManager() 
{
  factoryOn = false;

// Initialization
// histograms
  //for (G4int k=0; k<MaxHisto; k++) fHistId[k] = 0;

// Initialization ntuple
  for (G4int k=0; k<MaxNtCol; k++) fNtColId[k] = 0;

  //h10 = 0;
  //h20 = 0;
}

AnalysisManager::~AnalysisManager() 
{
}

void AnalysisManager::book() 
{ 
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  
  manager->SetVerboseLevel(1);
 
  // Create a root file
  G4String fileName = "BDUMP.root";

  // Create directories  
  //manager->SetNtupleDirectoryName("GeHP");

  G4bool fileOpen = manager->OpenFile(fileName);
  if (!fileOpen) {
    G4cout << "\n---> HistoManager::book(): cannot open " 
           << fileName << G4endl;
    return;
  }

  manager ->SetFirstH1Id(1);
  manager->SetFirstNtupleId(1);

  //Create Primary Energy Ntuple
  manager -> CreateNtuple("101", "Source Spectra");
  fNtColId[0] = manager -> CreateNtupleDColumn("Ek");
  manager -> FinishNtuple();

  
  manager -> CreateNtuple("102", "CR39-1");
  fNtColId[1] = manager -> CreateNtupleDColumn("edep");
  manager -> FinishNtuple();

  manager -> CreateNtuple("103", "CR39-2");
  fNtColId[1] = manager -> CreateNtupleDColumn("edep");
  manager -> FinishNtuple();

 
  //creating a ntuple, containing the information about secondary particles
  //manager -> CreateNtuple("103", "secondary");
  //fNtColId[2] = manager -> CreateNtupleDColumn("AA");
  //fNtColId[3] = manager -> CreateNtupleDColumn("ZZ");
  //fNtColId[4] = manager -> CreateNtupleDColumn("KE");
  //manager -> FinishNtuple();

  manager -> CreateH1("1", "Source Spectra", 1000, 2*keV, 2*MeV);
  manager -> CreateH1("2", "Detector Spectra", 1000, 2*keV, 2*MeV);
  factoryOn = true;    
}


void AnalysisManager::SetPrimaryEnergy(G4double energy)
{
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager -> FillH1(1, energy);
  manager -> FillNtupleDColumn(1, fNtColId[0], energy);
  manager -> AddNtupleRow(1); 
}

void AnalysisManager::StoreEnergyDeposition(G4double edep, int sd)
{
  if (edep > 1.1*MeV) {
    G4AnalysisManager* manager = G4AnalysisManager::Instance();
    manager -> FillNtupleDColumn(sd+1, fNtColId[1], edep);
    manager -> AddNtupleRow(2); 
    manager -> FillH1(sd, edep);
    G4cout << "Stored Edep " << edep << "Det: " << sd << G4endl;
  } else {
    G4cout << "Ign Edep " << edep << G4endl;
  }
  
  
}

void AnalysisManager::FillSecondaries(G4int AA, G4double charge, G4double energy)
{

//  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  G4cout << "Secondary " << AA << " charge " << charge << " energy " << energy << G4endl;
//  manager -> FillNtupleDColumn(3, fNtColId[2], AA);
//  manager -> FillNtupleDColumn(3, fNtColId[3], charge);
//  manager -> FillNtupleDColumn(3, fNtColId[4], energy);
 // manager -> AddNtupleRow(3);  
}
 
void AnalysisManager::finish() 
{   
 if (factoryOn) 
   {
    G4AnalysisManager* manager = G4AnalysisManager::Instance();    
    manager -> Write();
    manager -> CloseFile();  
      
    delete G4AnalysisManager::Instance();
    factoryOn = false;
   }
}












