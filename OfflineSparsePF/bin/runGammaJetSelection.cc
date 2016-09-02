//================================================================================================
//
// gamma+jets selection
//
// Input arguments:
//   argv[1] => input bacon file name
//   argv[2] => output bacon bits file name
//   argv[3] => dataset type: "mczsig", "mczbkg", "mcbkg", "data"
//   argv[4] => JSON file for run-lumi filtering of data, specify "none" for MC or no filtering
//   argv[5] => cross section (pb), ignored for data
//
//________________________________________________________________________________________________

// bacon object headers
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

// JSON file parser
#include "BaconAna/Utils/interface/RunLumiRangeMap.h"

// ROOT headers
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TMath.h>

// Other C++ headers
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>


//=== FUNCTION DECLARATIONS ======================================================================================

// Object selection
bool passPhoSel(const baconhep::TPhoton *photon, const double rho);
bool passJetSel(const baconhep::TJet *jet);

//=== MAIN =======================================================================================================

int main(int argc, char **argv)
{
  //--------------------------------------------------------------------------------------------------------------
  // Settings and constants
  //==============================================================================================================

  // handle input arguments
  const std::string infilename   = argv[1];
  const std::string outfilename  = argv[2];
  const std::string dstypename   = argv[3];
  const std::string jsonfilename = argv[4];
  const double      xsec         = atof(argv[5]);

  // determine dataset type
  enum {
    kMCZSig=1,  // Z->ee signal from MC
    kMCZBkg,    // Z->ll background from MC
    kMCBkg,     // other background from MC
    kData       // data
  };

  unsigned int dstype=0;
  if     (dstypename.compare("mcgammasig")==0) { dstype = kMCGammaSig; }
  else if(dstypename.compare("mcgammabkg")==0) { dstype = kMCGammaBkg; }
  else if(dstypename.compare("mcbkg")==0)  { dstype = kMCBkg; }
  else if(dstypename.compare("data")==0)   { dstype = kData; }
  assert(dstype>0);

  // Trigger bits mapping file
  std::string trigfilename = getenv("CMSSW_BASE");
  trigfilename += "/src/BaconAna/DataFormats/data/HLTFile_25ns";
  baconhep::TTrigger trigger(trigfilename);

  // Cuts
  const double PHO_PT_CUT  = 30;
  const double PHO_ETA_CUT = 2.1;

  //
  // constants
  //
  const int    PHO_PDGID   = 22;

  // Print summary of selection cuts
  std::cout << " ===== Cuts ===== " << std::endl;
  std::cout << " -- Photon: ";
  std::cout << " pT > "     << PHO_PT_CUT;
  std::cout << ", |eta| < " << PHO_ETA_CUT << std::endl;
  std::cout << std::endl;


  //--------------------------------------------------------------------------------------------------------------
  // Set up output file
  //==============================================================================================================

  unsigned int runNum, lumiSec, evtNum;       // event ID
  unsigned int npv;                           // number of PV
  unsigned int njets, njetsc;                 // jet multiplicity
  float        npu;                           // mean expected PU
  float        scale1fb;                      // event weight per 1/fb for xsec normalization
  float        lhew;                          // LHE weight
  float        pfmet, pfmetphi;               // MET
  float        pfmetraw, pfmetrawphi;
  float        trkmet, trkmetphi;
  float        puppet, puppetphi;

  TLorentzVector *pho=0;
  TLorentzVector *jet1=0, *jetc1=0;

  TFile *outFile = new TFile(outfilename.c_str(),"RECREATE");
  TH1D hTotalEvents("TotalEvents","TotalEvents",1,-10,10);
  TH1D hTotalWeights("TotalWeights","TotalWeights",1,-10,10);
  TTree *outTree = new TTree("Events","Events");
  
  outTree->Branch("runNum",      &runNum,      "runNum/i");
  outTree->Branch("lumiSec",     &lumiSec,     "lumiSec/i");
  outTree->Branch("evtNum",      &evtNum,      "evtNum/i");
  outTree->Branch("npv",         &npv,         "npv/i");
  outTree->Branch("njets",       &njets,       "njets/i");
  outTree->Branch("njetsc",      &njetsc,      "njetsc/i");
  outTree->Branch("npu",         &npu,         "npu/F");
  outTree->Branch("scale1fb",    &scale1fb,    "scale1fb/F");
  outTree->Branch("lhew",        &lhew,        "lhew/F");
  outTree->Branch("pfmet",       &pfmet,       "pfmet/F");
  outTree->Branch("pfmetphi",    &pfmetphi,    "pfmetphi/F");
  outTree->Branch("pfmetraw",    &pfmetraw,    "pfmetraw/F");
  outTree->Branch("pfmetrawphi", &pfmetrawphi, "pfmetrawphi/F");
  outTree->Branch("trkmet",      &trkmet,      "trkmet/F");
  outTree->Branch("trkmetphi",   &trkmetphi,   "trkmetphi/F");
  outTree->Branch("puppet",      &puppet,      "puppet/F");
  outTree->Branch("puppetphi",   &puppetphi,   "puppetphi/F");

  outTree->Branch("pho", "TLorentzVector", &pho);
  outTree->Branch("jet1",  "TLorentzVector", &jet1);
  outTree->Branch("jetc1", "TLorentzVector", &jetc1);


  //--------------------------------------------------------------------------------------------------------------
  // Process input file
  //==============================================================================================================
  baconhep::TEventInfo *info   = 0; TBranch *infoBr = 0;
  baconhep::TGenEventInfo *gen = 0; TBranch *genBr  = 0;
  TClonesArray *phoArr         = 0; TBranch *phoBr  = 0;
  TClonesArray *jetArr         = 0; TBranch *jetBr  = 0;
  TClonesArray *pvArr          = 0; TBranch *pvBr   = 0;
  TClonesArray *genParArr      = 0; TBranch *genParBr = 0;

  std::cout << "Processing " << infilename << "..." << std::endl;    
  TFile *infile    = TFile::Open(infilename.c_str()); assert(infile);
  TTree *eventTree = (TTree*)infile->Get("Events");   assert(eventTree);
  
  hTotalEvents.Add((TH1D*)infile->Get("TotalEvents"));
  if(dstype==kData) {
    hTotalWeights.Add((TH1D*)infile->Get("TotalEvents"));
  }

  eventTree->SetBranchAddress("Info",     &info,   &infoBr);
  eventTree->SetBranchAddress("Photon",   &phoArr, &phoBr);
  eventTree->SetBranchAddress("AK4CHS",   &jetArr, &jetBr);
  eventTree->SetBranchAddress("PV",       &pvArr,  &pvBr);
  if(dstype!=kData) {
    eventTree->SetBranchAddress("GenParticle", &genParArr, &genParBr);
    eventTree->SetBranchAddress("GenEvtInfo", &gen, &genBr);
    for(unsigned int ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      genBr->GetEntry(ientry);
      hTotalWeights.Fill(1, gen->weight);
    }
  }


  double weight = 1;  // event weight

  // Set up object to handle good run-lumi filtering if necessary
  bool hasJSON = false;
  RunLumiRangeMap rlrm;
  if(jsonfilename.compare("none")!=0) {
    rlrm.AddJSONFile(jsonfilename);
    hasJSON = true;
  }

  //
  // loop over events
  //
  for(unsigned int ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    infoBr->GetEntry(ientry);

    if(dstype==kData) {
      if(hasJSON) {
        // JSON filter
        RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);
        if(!rlrm.HasRunLumi(rl)) continue;
      }
    } else {
      if(dstype==kMCGammaSig || dstype==kMCGammaBkg) {
        genParArr->Clear();
        genParBr->GetEntry(ientry);

        bool isSig=false;
        for(int k=0; k<genParArr->GetEntriesFast(); k++) {
          const baconhep::TGenParticle *genp = (baconhep::TGenParticle*)genParArr->At(k);
          if(genp->pdgId==PHO_PDGID && genp->status==23) {
            isSig=true;
            break;
          }
          if(genp->pdgId==PHO_PDGID && genp->status==1) {
            int iparent = genp->parent;
            if(iparent>=0 && ((baconhep::TGenParticle*)genParArr->At(iparent))->pdgId==Z_PDGID) {
              isSig=true;
              break;
            }
          }
        }
        if(dstype==kMCGammaSig && !isSig) continue;
        if(dstype==kMCGammaBkg &&  isSig) continue;
      }

      genBr->GetEntry(ientry);
      if(xsec>0) {
        weight = 1000.*xsec*gen->weight/hTotalWeights.Integral();
      }
    }

    //
    // primary vertex requirement
    //
    if(!(info->hasGoodPV)) continue;

    // Read in PV branch to get number of primary vertices
    pvArr->Clear();
    pvBr->GetEntry(ientry);

    //
    // trigger requirement
    //
    std::string trigName = "HLT_Photon120_v*";
    if(!trigger.pass(trigName, info->triggerBits)) continue;

    //
    // photon selection: select isolated photon
    //                     passing kinematic cuts
    //
    phoArr->Clear();
    phoBr->GetEntry(ientry);

    for(int i=0; i<phoArr->GetEntriesFast(); i1++) {
      const baconhep::TPhoton *pho = (baconhep::TPhoton*)phoArr->At(i);
      
      if(pho->pt           <= PHO_PT_CUT) continue;
      if(fabs(pho->eta)    >= PHO_ETA_CUT) continue;
      if(!passPhoSel(pho, info->rhoIso))   continue;
      
      TLorentzVector vPho;
      vPho.SetPtEtaPhiM(pho->pt, pho->eta, pho->phi, 0);
            
      ///// We have a photon candidate! Hurray! /////
      
      // jet counting
      jetArr->Clear();
      jetBr->GetEntry(ientry);
      unsigned int nJets=0, nJetsC=0;
      TLorentzVector vJet1, vJetC1;
      for(int i=0; i<jetArr->GetEntriesFast(); i++) {
	const baconhep::TJet *jet = (baconhep::TJet*)jetArr->At(i);
	
	TLorentzVector vJet;
	vJet.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
	
	if(vJet.Pt()        <= 30)  continue;
	if(fabs(vJet.Eta()) >= 4.0) continue;
	if(!passJetSel(jet))        continue;
	
	nJets++;
	if(vJet.Pt() > vJet1.Pt()) { vJet1 = vJet; }
	
	if(fabs(vJet.Eta()) < 2.4) {
	  nJetsC++;
	  if(vJet.Pt() > vJetC1.Pt()) { vJetC1 = vJet; }
	}
      }
      
      //
      // Fill output tree
      //
      runNum      = info->runNum;
      lumiSec     = info->lumiSec;
      evtNum      = info->evtNum;
      npv         = pvArr->GetEntriesFast();
      njets       = nJets;
      njetsc      = nJetsC;
      npu         = info->nPUmean;
      scale1fb    = weight;
      lhew        = (dstype!=kData) ? gen->weight : 1.;
      pfmet       = info->pfMETC;
      pfmetphi    = info->pfMETCphi;
      pfmetraw    = info->pfMET;
      pfmetrawphi = info->pfMETphi;
      trkmet      = info->trkMET;
      trkmetphi   = info->trkMETphi;
      puppet      = info->puppET;
      puppetphi   = info->puppETphi;
      pho         = &vPho;
      
      // leading jet info
      jet1  = &vJet1;
      jetc1 = &vJetC1;
      
      outTree->Fill();
      
    }
  }
  
  delete infile;
  infile=0, eventTree=0;
  
  outFile->Write();
  outFile->Close();
  
  delete info;
  delete gen;
  delete phoArr;
  delete jetArr;
  delete pvArr;
  delete genParArr;

  return 0;
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
bool passPhoSel(const baconhep::TPhoton *photon, const double rho)
{

  // Tight cut-based photon selection
  // (https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2#SPRING15_selections_25_ns)
  if(!(photon->passElectronVeto)) return false;  // conversion safe electron veto


  double chHadIso  = TMath::Max(photon->chHadIso  - rho*phoEffArea(photon->scEta, 0), (double)0.);
  double neuHadIso = TMath::Max(photon->neuHadIso - rho*phoEffArea(photon->scEta, 1), (double)0.);
  double phoIso    = TMath::Max(photon->gammaIso  - rho*phoEffArea(photon->scEta, 2), (double)0.);

  if(fabs(photon->scEta) <= 1.479) {
    if(photon->sthovere > 0.05)                                                           return false;
    if(photon->sieie    > 0.0100)                                                         return false;
    if(photon->sieie    > 0.0102)                                                         return false;
    if(chHadIso         > 0.76)                                                           return false;
    if(neuHadIso        > 0.97 + 0.014*(photon->pt) + 0.000019*(photon->pt)*(photon->pt)) return false;
    if(phoIso           > 0.08 + 0.0053*photon->pt)                    return false;

  } else {
    if(photon->sthovere > 0.05)                                                            return false;
    if(photon->sieie    > 0.0268)                                                          return false;
    if(chHadIso         > 0.56)                                                            return false;
    if(neuHadIso        > 2.09 + 0.0139*(photon->pt) + 0.000025*(photon->pt)*(photon->pt)) return false;
    if(phoIso           > 0.83 + 0.0034*photon->pt)                                        return false;
  }

  return true;
}


//--------------------------------------------------------------------------------------------------
bool passJetSel(const baconhep::TJet *jet)
{
  // Loose PFJet ID
  // see twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
  if(fabs(jet->eta)<=3.0) {
    if(jet->neuHadFrac >= 0.99) return false;
    if(jet->neuEmFrac  >= 0.99) return false;
    if(jet->nParticles <= 1)    return false;
    //if(jet->muonFrac   >= 0.80) return false;
    if(fabs(jet->eta)<=2.4) {
      if(jet->chHadFrac == 0)    return false;
      if(jet->nCharged  == 0)    return false;
      if(jet->chEmFrac  >= 0.99) return false;
    }
  } else {
    if(jet->neuEmFrac >= 0.90) return false;
    if(jet->nNeutrals <= 10)   return false;
  }
/*
  // PU Jet ID
  if     (0    <= fabs(jet->eta) && fabs(jet->eta) < 2.5  && jet->mva < -0.63) return false;
  else if(2.5  <= fabs(jet->eta) && fabs(jet->eta) < 2.75 && jet->mva < -0.60) return false;
  else if(2.75 <= fabs(jet->eta) && fabs(jet->eta) < 3    && jet->mva < -0.55) return false;
  else if(3    <= fabs(jet->eta) && fabs(jet->eta) < 5    && jet->mva < -0.45) return false;
*/
  return true;
}
