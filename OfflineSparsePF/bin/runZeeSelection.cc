//================================================================================================
//
// Z->ee selection
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
bool passEleSel(const baconhep::TElectron *electron, const double rho);
double eleEffArea(const double eta);
bool passJetSel(const baconhep::TJet *jet);
/*
// energy scale corrections for data
double getScaleCorr(const double eta)
{
  if     (fabs(eta) < 0.4) { return 1.00387; }
  else if(fabs(eta) < 0.8) { return 1.00158; }
  else if(fabs(eta) < 1.2) { return 1.00181; }
  else if(fabs(eta) < 1.4) { return 0.997462; }
  else if(fabs(eta) < 1.6) { return 1.00517; }
  else if(fabs(eta) < 2.1) { return 1.00645; }
  else                     { return 1.02714; }
}
*/

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
  if     (dstypename.compare("mczsig")==0) { dstype = kMCZSig; }
  else if(dstypename.compare("mczbkg")==0) { dstype = kMCZBkg; }
  else if(dstypename.compare("mcbkg")==0)  { dstype = kMCBkg; }
  else if(dstypename.compare("data")==0)   { dstype = kData; }
  assert(dstype>0);

  // Trigger bits mapping file
  std::string trigfilename = getenv("CMSSW_BASE");
  trigfilename += "/src/BaconAna/DataFormats/data/HLTFile_25ns";
  baconhep::TTrigger trigger(trigfilename);

  // Cuts
  const double ELE_PT_CUT  = 30;
  const double ELE_ETA_CUT = 2.1;
  const double ZMASSLOW    = 40;
  const double ZMASSHIGH   = 200;

  //
  // constants
  //
  const double ELE_MASS  = 0.000511;
  const int    ELE_PDGID = 11;
  const int    Z_PDGID   = 23;

  // Print summary of selection cuts
  std::cout << " ===== Cuts ===== " << std::endl;
  std::cout << " -- Electron: ";
  std::cout << " pT > "     << ELE_PT_CUT;
  std::cout << ", |eta| < " << ELE_ETA_CUT << std::endl;
  std::cout << " -- Z mass window: " << ZMASSLOW << " < M < " << ZMASSHIGH << std::endl;
  std::cout << std::endl;


  //--------------------------------------------------------------------------------------------------------------
  // Set up output file
  //==============================================================================================================

  unsigned int runNum, lumiSec, evtNum;       // event ID
  unsigned int npv;                           // number of PV
  unsigned int njets, njetsc;                 // jet multiplicity
  int          q1, q2;                        // charge
  float        npu;                           // mean expected PU
  float        scale1fb;                      // event weight per 1/fb for xsec normalization
  float        lhew;                          // LHE weight
  float        pfmet, pfmetphi;               // MET
  float        pfmetraw, pfmetrawphi;
  float        trkmet, trkmetphi;
  float        puppet, puppetphi;
  TLorentzVector *lep1=0, *lep2=0, *dilep=0;  // 4-vectors
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

  outTree->Branch("q1", &q1, "q1/I");
  outTree->Branch("lep1", "TLorentzVector", &lep1);

  outTree->Branch("q2", &q2, "q2/I");
  outTree->Branch("lep2", "TLorentzVector", &lep2);

  outTree->Branch("dilep", "TLorentzVector", &dilep);

  outTree->Branch("jet1",  "TLorentzVector", &jet1);
  outTree->Branch("jetc1", "TLorentzVector", &jetc1);


  //--------------------------------------------------------------------------------------------------------------
  // Process input file
  //==============================================================================================================
  baconhep::TEventInfo *info   = 0; TBranch *infoBr = 0;
  baconhep::TGenEventInfo *gen = 0; TBranch *genBr  = 0;
  TClonesArray *eleArr         = 0; TBranch *eleBr  = 0;
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
  eventTree->SetBranchAddress("Electron", &eleArr, &eleBr);
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
      if(dstype==kMCZSig || dstype==kMCZBkg) {
        genParArr->Clear();
        genParBr->GetEntry(ientry);

        bool isSig=false;
        for(int k=0; k<genParArr->GetEntriesFast(); k++) {
          const baconhep::TGenParticle *genp = (baconhep::TGenParticle*)genParArr->At(k);
          if(genp->pdgId==ELE_PDGID && genp->status==23) {
            isSig=true;
            break;
          }
          if(genp->pdgId==ELE_PDGID && genp->status==1) {
            int iparent = genp->parent;
            if(iparent>=0 && ((baconhep::TGenParticle*)genParArr->At(iparent))->pdgId==Z_PDGID) {
              isSig=true;
              break;
            }
          }
        }
        if(dstype==kMCZSig && !isSig) continue;
        if(dstype==kMCZBkg &&  isSig) continue;
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
    std::string trigName = "HLT_Ele27_eta2p1_WP75_Gsf_v*";
    if(dstype==kData) { trigName = "HLT_Ele27_eta2p1_WPLoose_Gsf_v*"; }
    if(!trigger.pass(trigName, info->triggerBits)) continue;

    //
    // electron selection: select pairs of oppositely charged dielectrons
    //                     passing kinematic cuts and mass window
    //
    eleArr->Clear();
    eleBr->GetEntry(ientry);

    for(int i1=0; i1<eleArr->GetEntriesFast(); i1++) {
      const baconhep::TElectron *ele1 = (baconhep::TElectron*)eleArr->At(i1);

      double scalecorr1 = 1;//(dstype==kData) ? getScaleCorr(ele1->eta) : 1;
      if(ele1->pt*scalecorr1 <= ELE_PT_CUT) continue;
      if(fabs(ele1->eta)    >= ELE_ETA_CUT) continue;
      if(!passEleSel(ele1, info->rhoIso))   continue;

      TLorentzVector vEle1;
      vEle1.SetPtEtaPhiM(ele1->pt*scalecorr1, ele1->eta, ele1->phi, ELE_MASS);

      for(int i2=i1+1; i2<eleArr->GetEntriesFast(); i2++) {
        const baconhep::TElectron *ele2 = (baconhep::TElectron*)eleArr->At(i2); 
        if(ele1->q == ele2->q) continue;

        double scalecorr2 = 1;//(dstype==kData) ? getScaleCorr(ele2->eta) : 1;
        if(ele2->pt*scalecorr2 <= ELE_PT_CUT) continue;
        if(fabs(ele2->eta)    >= ELE_ETA_CUT) continue;
        if(!passEleSel(ele2, info->rhoIso))   continue;

        TLorentzVector vEle2;
        vEle2.SetPtEtaPhiM(ele2->pt*scalecorr2, ele2->eta, ele2->phi, ELE_MASS);

        // mass window
        TLorentzVector vDilep = vEle1 + vEle2;
        if(vDilep.M()<ZMASSLOW || vDilep.M()>ZMASSHIGH) continue;

        // trigger object matching
        if(!trigger.passObj(trigName, 1, ele1->hltMatchBits) &&
           !trigger.passObj(trigName, 1, ele2->hltMatchBits))
          continue;

        ///// We have a Z candidate! Hurray! /////

        // jet counting
        jetArr->Clear();
        jetBr->GetEntry(ientry);
        unsigned int nJets=0, nJetsC=0;
        TLorentzVector vJet1, vJetC1;
        for(int i=0; i<jetArr->GetEntriesFast(); i++) {
          const baconhep::TJet *jet = (baconhep::TJet*)jetArr->At(i);

          TLorentzVector vJet;
          vJet.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);

          if(vJet.DeltaR(vEle1) < 0.4) continue;
          if(vJet.DeltaR(vEle2) < 0.4) continue;

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
        dilep       = &vDilep;
  
        // leading jet info
        jet1  = &vJet1;
        jetc1 = &vJetC1;

        // leading electron info
        q1   = (vEle1.Pt()>vEle2.Pt()) ? ele1->q : ele2->q;
        lep1 = (vEle1.Pt()>vEle2.Pt()) ? &vEle1 : &vEle2;

        // trailing electron info
        q2   = (vEle1.Pt()>vEle2.Pt()) ? ele2->q : ele1->q;
        lep2 = (vEle1.Pt()>vEle2.Pt()) ? &vEle2 : &vEle1;

        outTree->Fill();
      }
    }
  }
  
  delete infile;
  infile=0, eventTree=0;
  
  outFile->Write();
  outFile->Close();
  
  delete info;
  delete gen;
  delete eleArr;
  delete jetArr;
  delete pvArr;
  delete genParArr;

  return 0;
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
bool passEleSel(const baconhep::TElectron *electron, const double rho)
{
  // Spring15 25ns V1 Tight cut-based selection
  // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns

  if(electron->isConv) return false;

  double iso = electron->chHadIso + TMath::Max( 0.0,(electron->gammaIso + electron->neuHadIso - rho*eleEffArea(electron->eta)) );

  if(fabs(electron->scEta)<1.479) {
    if(iso >= 0.0354*(electron->pt)) return false;

    if(electron->sieie              >= 0.0101)                       return false;
    if(fabs(electron->dEtaIn)       >= 0.00926)                      return false;
    if(fabs(electron->dPhiIn)       >= 0.0336)                       return false;
    if(electron->hovere             >= 0.0597)                       return false;
    if(fabs(1.0 - electron->eoverp) >= 0.012*(electron->ecalEnergy)) return false;
    if(fabs(electron->d0)           >= 0.0111)                       return false;
    if(fabs(electron->dz)           >= 0.0466)                       return false;
    if(electron->nMissingHits       >  2)                            return false;

  } else {
    if(iso >= 0.0646*(electron->pt)) return false;

    if(electron->sieie              >= 0.0279)                         return false;
    if(fabs(electron->dEtaIn)       >= 0.00724)                        return false;
    if(fabs(electron->dPhiIn)       >= 0.0918)                         return false;
    if(electron->hovere             >= 0.0615)                         return false;
    if(fabs(1.0 - electron->eoverp) >= 0.00999*(electron->ecalEnergy)) return false;
    if(fabs(electron->d0)           >= 0.0351)                         return false;
    if(fabs(electron->dz)           >= 0.417)                          return false;
    if(electron->nMissingHits       >  1)                              return false;
  }

  return true;
}

//--------------------------------------------------------------------------------------------------
double eleEffArea(const double eta)
{
  // effective area for PU correction
  // (see slide 12 of https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf)

  if     (fabs(eta) >= 0.0   && fabs(eta) < 1.0)   { return 0.1752; }
  else if(fabs(eta) >= 1.0   && fabs(eta) < 1.479) { return 0.1862; }
  else if(fabs(eta) >= 1.479 && fabs(eta) < 2.0)   { return 0.1411; }
  else if(fabs(eta) >= 2.0   && fabs(eta) < 2.2)   { return 0.1534; }
  else if(fabs(eta) >= 2.2   && fabs(eta) < 2.3)   { return 0.1903; }
  else if(fabs(eta) >= 2.3   && fabs(eta) < 2.4)   { return 0.2243; }
  else                                             { return 0.2687; }
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
