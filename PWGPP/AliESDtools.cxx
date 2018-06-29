/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
    author marian.ivanov@cern.ch
*/


/*
  TFile * f = TFile::Open("/lustre/nyx/alice/users/marsland/alice-tpc-notes/JIRA/ATO-436/data/LHC15o_pass1_esds/alice/data/2015/LHC15o/000246272/pass1/15000246272021.7701/AliESDs.root");
  TTree * tree = (TTree*)f->Get("esdTree");
  // tree = AliXRDPROOFtoolkit::MakeChainRandom("esd.list","esdTree",0,10)
  .L $AliPhysics_SRC/PWGPP/AliESDtools.cxx+
  AliESDtools tools;
  tools.Init(tree);
   //tree->GetEntry(3);
  //tree->Scan("Tracks@.GetEntries():SPDVertex.fNContributors:Entry$","SPDVertex.fNContributors>100&&Tracks@.GetEntries()/PrimaryVertex.fNContributors>10")
  tools->CacheTPCEventInformation()

  /// trigger exceptional pileup
  tree->Draw(">>entryList","SPDVertex.fNContributors>100&&Tracks@.GetEntries()/PrimaryVertex.fNContributors>10","entrylist");
   tree->SetEntryList(entryList)
   tree->Scan("AliESDtools::GetTrackMatchEff(0,0):AliESDtools::GetTrackCounters(0,0):AliESDtools::GetTrackCounters(4,0):AliESDtools::GetMeanHisTPCVertexA():AliESDtools::GetMeanHisTPCVertexC():Entry$","AliESDtools::SCalculateEventVariables(Entry$)")
*/


#include "TStopwatch.h"
#include "TTree.h" 
#include "TChain.h"
#include "TVectorF.h"
#include "AliStack.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TTreeStream.h"
#include "AliRunLoader.h"
#include "AliTrackReference.h"
#include "AliExternalTrackParam.h"
#include "AliHelix.h"
#include "TCut.h"
#include "AliTreePlayer.h"
#include "THn.h"
#include "TF3.h"
#include "TStatToolkit.h"
#include <stdarg.h>
#include "AliNDLocalRegression.h"
#include "AliESDEvent.h"
#include "AliESDtools.h"
#include "AliTrackerBase.h"

ClassImp(AliESDtools)
AliESDtools*  AliESDtools::fgInstance;

AliESDtools::AliESDtools():
  fVerbose(0),
  fESDtree(NULL),
  fEvent(NULL),
  fHisTPCVertexA(nullptr),
  fHisTPCVertexC(nullptr),
  fHisTPCVertexACut(nullptr),
  fHisTPCVertexCCut(nullptr),
  fHisTPCVertex(nullptr),
  fHistPhiTPCcounterA(nullptr),         // helper histogram phi counteres
  fHistPhiTPCcounterC(nullptr),         // helper histogram for TIdentity tree
  fHistPhiTPCcounterAITS(nullptr),      // helper histogram for TIdentity tree
  fHistPhiTPCcounterCITS(nullptr),      // helper histogram for TIdentity tree
  fHistPhiITScounterA(nullptr),         // helper histogram for TIdentity tree
  fHistPhiITScounterC(nullptr),         // helper histogram for TIdentity tree
  fCacheTrackCounters(nullptr),         // track counter
  fCacheTrackTPCCountersZ(nullptr),     // track counter
  fCacheTrackdEdxRatio(nullptr),        // dedx info counter
  fCacheTrackNcl(nullptr),              // ncl counter
  fCacheTrackChi2(nullptr),             // chi2 counter
  fCacheTrackMatchEff(nullptr),         // matchEff counter
  fLumiGraph(nullptr),                  // graph for the interaction rate info for a run
  fGlobalID(0),
  fStreamer(nullptr)
{
  fgInstance=this;
}

AliESDtools::AliESDtools(const AliESDtools &tools) :
  TNamed(tools),
  fVerbose(tools.fVerbose),
  fESDtree(tools.fESDtree),
  fEvent(tools.fEvent),
  fHisTPCVertexA(tools.fHisTPCVertexA),
  fHisTPCVertexC(tools.fHisTPCVertexC),
  fHisTPCVertexACut(tools.fHisTPCVertexACut),
  fHisTPCVertexCCut(tools.fHisTPCVertexCCut),
  fHisTPCVertex(tools.fHisTPCVertex),
  fHistPhiTPCcounterA(tools.fHistPhiTPCcounterA),         // helper histogram phi counteres
  fHistPhiTPCcounterC(tools.fHistPhiTPCcounterC),         // helper histogram for TIdentity tree
  fHistPhiTPCcounterAITS(tools.fHistPhiTPCcounterAITS),   // helper histogram for TIdentity tree
  fHistPhiTPCcounterCITS(tools.fHistPhiTPCcounterCITS),   // helper histogram for TIdentity tree
  fHistPhiITScounterA(tools.fHistPhiITScounterA),         // helper histogram for TIdentity tree
  fHistPhiITScounterC(tools.fHistPhiITScounterA),         // helper histogram for TIdentity tree
  fCacheTrackCounters(tools.fCacheTrackCounters),         // track counter
  fCacheTrackTPCCountersZ(tools.fCacheTrackTPCCountersZ), // track counter
  fCacheTrackdEdxRatio(tools.fCacheTrackdEdxRatio),       // dedx info counter
  fCacheTrackNcl(tools.fCacheTrackNcl),                   // ncl counter
  fCacheTrackChi2(tools.fCacheTrackChi2),                 // chi2 counter
  fCacheTrackMatchEff(tools.fCacheTrackMatchEff),         // matchEff counter
  fLumiGraph(tools.fLumiGraph),                           // graph for the interaction rate info for a run
  fGlobalID(tools.fGlobalID),
  fStreamer(tools.fStreamer)
{
  // copy constructor
}

AliESDtools& AliESDtools::operator=(const AliESDtools &rhs)
{
  // assignment operator
  if (this != &rhs) {
    TNamed::operator=(rhs);
    fVerbose                = rhs.fVerbose; 
    fESDtree                = rhs.fESDtree; 
    fEvent                  = rhs.fEvent; 
    fHisTPCVertexA          = rhs.fHisTPCVertexA; 
    fHisTPCVertexC          = rhs.fHisTPCVertexC; 
    fHisTPCVertexACut       = rhs.fHisTPCVertexACut; 
    fHisTPCVertexCCut       = rhs.fHisTPCVertexCCut; 
    fHisTPCVertex           = rhs.fHisTPCVertex; 
    fHistPhiTPCcounterA     = rhs.fHistPhiTPCcounterA; 
    fHistPhiTPCcounterC     = rhs.fHistPhiTPCcounterC; 
    fHistPhiTPCcounterAITS  = rhs.fHistPhiTPCcounterAITS; 
    fHistPhiTPCcounterCITS  = rhs.fHistPhiTPCcounterCITS; 
    fHistPhiITScounterA     = rhs.fHistPhiITScounterA; 
    fHistPhiITScounterC     = rhs.fHistPhiITScounterC; 
    fCacheTrackCounters     = rhs.fCacheTrackCounters; 
    fCacheTrackTPCCountersZ = rhs.fCacheTrackTPCCountersZ; 
    fCacheTrackdEdxRatio    = rhs.fCacheTrackdEdxRatio; 
    fCacheTrackNcl          = rhs.fCacheTrackNcl; 
    fCacheTrackChi2         = rhs.fCacheTrackChi2; 
    fCacheTrackMatchEff     = rhs.fCacheTrackMatchEff; 
    fLumiGraph              = rhs.fLumiGraph; 
    fGlobalID               = rhs.fGlobalID;
    fStreamer               = rhs.fStreamer;
  }
  return *this;
}


void AliESDtools::Init(TTree *tree) {
  AliESDtools & tools = *this;
  if (tools.fESDtree) delete tools.fESDtree;
  if (!tools.fEvent) tools.fEvent = new AliESDEvent();
  tools.fESDtree = tree;
  tools.fEvent->ReadFromTree(tree);
  if (fHisTPCVertexA == nullptr) {
    tools.fHisTPCVertexA = new TH1F("hisTPCZA", "hisTPCZA", 1000, -250, 250);
    tools.fHisTPCVertexC = new TH1F("hisTPCZC", "hisTPCZC", 1000, -250, 250);
    tools.fHisTPCVertex = new TH1F("hisTPCZ", "hisTPCZ", 1000, -250, 250);
    tools.fHisTPCVertexACut = new TH1F("hisTPCZACut", "hisTPCZACut", 1000, -250, 250);
    tools.fHisTPCVertexCCut = new TH1F("hisTPCZCCut", "hisTPCZCCut", 1000, -250, 250);
    tools.fHisTPCVertex->SetLineColor(1);
    tools.fHisTPCVertexA->SetLineColor(2);
    tools.fHisTPCVertexC->SetLineColor(4);
    tools.fHisTPCVertexACut->SetLineColor(3);
    tools.fHisTPCVertexCCut->SetLineColor(6);
    tools.fCacheTrackCounters = new TVectorF(20);
    tools.fCacheTrackTPCCountersZ = new TVectorF(8);
    tools.fCacheTrackdEdxRatio = new TVectorF(20);
    tools.fCacheTrackNcl = new TVectorF(20);
    tools.fCacheTrackChi2 = new TVectorF(20);
    tools.fCacheTrackMatchEff = new TVectorF(20);
    // **************************** Event histograms **************************
    fHistPhiTPCcounterA = new TH1F("hPhiTPCcounterC", "control histogram to count tracks on the A side in phi ", 36, 0., 18.);
    fHistPhiTPCcounterC = new TH1F("hPhiTPCcounterA", "control histogram to count tracks on the C side in phi ", 36, 0., 18.);
    fHistPhiTPCcounterAITS = new TH1F("hPhiTPCcounterAITS", "control histogram to count tracks on the A side in phi ", 36, 0., 18.);
    fHistPhiTPCcounterCITS = new TH1F("hPhiTPCcounterCITS", "control histogram to count tracks on the C side in phi ", 36, 0., 18.);
    fHistPhiITScounterA = new TH1F("hPhiITScounterA", "control histogram to count tracks on the A side in phi ", 36, 0., 18.);
    fHistPhiITScounterC = new TH1F("hPhiITScounterC", "control histogram to count tracks on the C side in phi ", 36, 0., 18.);
  }
  fStreamer = new TTreeSRedirector("esdToolsDebug.root", "recreate");
}

void AliESDtools::DeInit() {
  if (fStreamer) {
    delete fStreamer;
  }
}

/// cache TPC event information
/// \return
Int_t AliESDtools::CacheTPCEventInformation(){
  AliESDtools &tools=*this;
  const Int_t kNCRCut=80;
  const Double_t kDCACut=5;
  const Float_t knTrackletCut=1.5;
  // FILL DCA histograms
  tools.fHisTPCVertexA->Reset();
  tools.fHisTPCVertexC->Reset();
  tools.fHisTPCVertexACut->Reset();
  tools.fHisTPCVertexCCut->Reset();
  tools.fHisTPCVertex->Reset();
  Int_t nTracks=tools.fEvent->GetNumberOfTracks();
  Int_t selected=0;
  for (Int_t iTrack=0; iTrack<nTracks; iTrack++){
    AliESDtrack * pTrack = tools.fEvent->GetTrack(iTrack);
    Float_t dcaxy,dcaz;
    if (pTrack== nullptr) continue;
    if (pTrack->IsOn(AliVTrack::kTPCin)==0) continue;
    if (pTrack->GetTPCClusterInfo(3,1)<kNCRCut) continue;
    pTrack->GetImpactParameters(dcaxy,dcaz);
    if (TMath::Abs(dcaxy)>kDCACut) continue;
    pTrack->SetESDEvent(fEvent);
    selected++;
    if ((pTrack->GetNumberOfTRDClusters()/20.+pTrack->GetNumberOfITSClusters())>knTrackletCut){
      tools.fHisTPCVertex->Fill(pTrack->GetTPCInnerParam()->GetZ());
      if (pTrack->GetTgl()>0) tools.fHisTPCVertexACut->Fill(pTrack->GetTPCInnerParam()->GetZ());
      if (pTrack->GetTgl()<0) tools.fHisTPCVertexCCut->Fill(pTrack->GetTPCInnerParam()->GetZ());
    }else{
      if (pTrack->GetTgl()>0) tools.fHisTPCVertexA->Fill(pTrack->GetTPCInnerParam()->GetZ());
      if (pTrack->GetTgl()<0) tools.fHisTPCVertexC->Fill(pTrack->GetTPCInnerParam()->GetZ());

    }
  }
  /*
  TPCVertexFit(tools.fHisTPCVertex);
  TPCVertexFit(tools.fHisTPCVertexA);
  TPCVertexFit(tools.fHisTPCVertexC);
  TPCVertexFit(tools.fHisTPCVertexACut);
  TPCVertexFit(tools.fHisTPCVertexCCut);
   */
  if (fVerbose&0x10) printf("%d\n",selected);//acheTPCEventInformation()
  //
  return selected;
}

void AliESDtools::TPCVertexFit(TH1F *hisVertex){
  //0.5 cm 1 bin
  // hisVertex=tools.fHisTPCVertexACut;
  TAxis * axis = hisVertex->GetXaxis();
  for (Int_t iBin=5; iBin<axis->GetNbins()-5; iBin++){
    Double_t median10=TMath::Median(10,&(hisVertex->GetArray()[iBin-5]));
    Double_t rms10=TMath::RMS(10,&(hisVertex->GetArray()[iBin-5]));
    Double_t val0=TMath::Mean(3,&(hisVertex->GetArray()[iBin-2]));
    Double_t val1=TMath::Mean(3,&(hisVertex->GetArray()[iBin-1]));
    Double_t val2=TMath::Mean(3,&(hisVertex->GetArray()[iBin+0]));
    if (val1>=val0 && val1>=val2 && val1>3+1.5*median10){
      Double_t xbin=axis->GetBinCenter(iBin);
      printf("Ibin %d\t%f\t%f\t%f\t%f\n", iBin,xbin,val1,median10,rms10);
      hisVertex->Fit("gaus","qnrsame+","qnr",xbin-2,xbin+2);
    }
  }
}

//
//
Int_t AliESDtools::CalculateEventVariables() {
  //AliVEvent *event=InputEvent();
  CacheTPCEventInformation();
  fCacheTrackCounters->Zero();   // track counter
  fCacheTrackCounters->Zero();   // track counter
  fCacheTrackdEdxRatio->Zero(); // dedx info counter
  fCacheTrackNcl->Zero();;       // ncl counter
  fCacheTrackChi2->Zero();;      // chi2 counter
  fCacheTrackMatchEff->Zero();  // matchEff counter
  //
  //
  if (fHistPhiTPCcounterA) fHistPhiTPCcounterA->Reset();
  if (fHistPhiTPCcounterC) fHistPhiTPCcounterC->Reset();
  if (fHistPhiTPCcounterAITS) fHistPhiTPCcounterAITS->Reset();
  if (fHistPhiTPCcounterCITS) fHistPhiTPCcounterCITS->Reset();
  if (fHistPhiITScounterA) fHistPhiITScounterA->Reset();
  if (fHistPhiITScounterC) fHistPhiITScounterC->Reset();
  //
  //
  const Int_t kNclTPCcut=60;
  const Float_t kTglCut=1.5;
  const Int_t kDCACut=5;  // 5 cm primary cut
  const Int_t kMindEdxClustersRegion=15;
  const Float_t kPtCut=0.100;
  //
  //
  ULong64_t orbitID = (ULong64_t) fEvent->GetOrbitNumber();
  ULong64_t bunchCrossID = (ULong64_t) fEvent->GetBunchCrossNumber();
  ULong64_t periodID = (ULong64_t) fEvent->GetPeriodNumber();
  fGlobalID = ( (periodID << 36ULL) || (orbitID << 12ULL) | (bunchCrossID) );
  //printf("orbit=%llu, bunchCross=%llu, period=%llu, gid=%llu\n", orbitID, bunchCrossID, periodID, fGlobalID);
  //
  //
  AliTPCdEdxInfo tpcdEdxInfo;
  for (Int_t itrack=0;itrack<fEvent->GetNumberOfTracks();++itrack)
  {   // Track loop

    AliESDtrack *track = fEvent->GetTrack(itrack);
    // AliESDPid *pid = track->GetDetPid();
    Double_t eta=-100., phiTPC=0.,sectorNumber=0.;
    Double_t tgl     = track->Pz()/track->Pt();
    Double_t phi  = track->Phi()-TMath::Pi();
    phi = track->GetParameterAtRadius(85,5,7);
    Double_t sectorNumbertmp = (9*phi/TMath::Pi()+18*(phi<0));
    if (track == NULL) continue;
    eta = track->Eta();
    if (TMath::Abs(eta)>0.9) continue;
    Bool_t isOnITS = track->IsOn(AliESDtrack::kITSrefit);
    Bool_t isOnTRD = track->IsOn(AliESDtrack::kTRDrefit);
    Bool_t isOnTPC = track->IsOn(AliESDtrack::kTPCrefit);
    //
    //
    if (track->GetInnerParam()) {
      phiTPC = track->GetInnerParam()->GetParameterAtRadius(85,5,7);
      sectorNumber = (9*phiTPC/TMath::Pi()+18*(phiTPC<0));
    }
    //
    // Count only ITS tracks
    if ( isOnITS && !isOnTPC ) {
      if (TMath::Abs(phi)>1e-10){
        if (tgl>0) fHistPhiITScounterA->Fill(sectorNumbertmp);
        if (tgl<0) fHistPhiITScounterC->Fill(sectorNumbertmp);
      }
    }
    //
    if (!track->GetInnerParam()) continue;
    if (track->IsOn(AliVTrack::kTPCout)==kFALSE)  continue;
    (*fCacheTrackCounters)[4]++;      // all TPC track with out flag
    // TPC track counters with DCAZ
    for (Int_t izCut=1; izCut<4; izCut++){
      Float_t impactParam[2];
      track->GetImpactParameters(impactParam[0],impactParam[1]);
      if (TMath::Abs(impactParam[0])>kDCACut) continue;
      if (TMath::Abs(track->GetInnerParam()->GetParameter()[1])<10.*(izCut+1.)) (*fCacheTrackTPCCountersZ)[izCut]++;
      if (TMath::Abs(impactParam[1])<10.*(izCut+1.)) (*fCacheTrackTPCCountersZ)[izCut+4]++;
    }
    //
    //
    Float_t dcaRPhi, dcaZ;
    track->GetImpactParameters(dcaRPhi, dcaZ);
    Int_t nclTPC    = track->GetTPCncls();
    Int_t nclITS    = track->GetITSNcls();
    Int_t nclTRD    = track->GetTRDncls();
    Int_t nclTOF    = track->IsOn(AliVTrack::kTOFout);
    if (nclTRD<1) nclTRD=-1;
    if (nclITS<1) nclITS=-1;
    //if (fNcl<1) fNcl=-1;
    Double_t chi2TPC = TMath::Sqrt(TMath::Abs(track->GetTPCchi2()/nclTPC));
    Double_t chi2ITS = TMath::Sqrt(TMath::Abs(track->GetITSchi2()));
    Double_t chi2TRD = TMath::Sqrt(TMath::Abs(track->GetTRDchi2()));
    Double_t ptot0   = track->GetP();
    Double_t qP      = track->Charge()/track->P();
    //
    //
    if (nclTPC<kNclTPCcut) continue;
    if (TMath::Abs(tgl)>kTglCut) continue;
    if (track->Pt()<kPtCut) continue;
    if (TMath::Abs(dcaRPhi)>kDCACut || TMath::Abs(dcaZ)>kDCACut) continue;
    (*fCacheTrackCounters)[5]++;
    if (TMath::Abs(phiTPC)>1e-10){
      if (tgl>0) fHistPhiTPCcounterA->Fill(sectorNumber);
      if (tgl<0) fHistPhiTPCcounterC->Fill(sectorNumber);
      if(isOnITS){
        if (tgl>0) fHistPhiTPCcounterAITS->Fill(sectorNumber);
        if (tgl<0) fHistPhiTPCcounterCITS->Fill(sectorNumber);
      }
    }
    //
    //
    Bool_t pileUpCut=  ( (nclITS>2) || (nclTRD>40));
    if (pileUpCut==kFALSE) continue;
    if (TMath::Min(chi2TPC,100.)<0) continue;
    (*fCacheTrackCounters)[1]++;

    //
    Bool_t itsOK=track->IsOn(AliVTrack::kITSout) && nclITS>2  && chi2ITS>0;
    Bool_t trdOK=track->IsOn(AliVTrack::kTRDout) && nclTRD>35 && chi2TRD>0;
    Bool_t tofOK=track->IsOn(AliVTrack::kTOFout);
    //
    (*fCacheTrackNcl)[4]+=track->GetTPCncls(0, 63);
    (*fCacheTrackNcl)[5]+=track->GetTPCncls(64, 127);
    (*fCacheTrackNcl)[6]+=track->GetTPCncls(128, 159);
    (*fCacheTrackNcl)[1] += nclTPC;
    (*fCacheTrackChi2)[1]+= (chi2TPC>0) ? TMath::Sqrt(chi2TPC):2;   // sometimes negative chi2?

    if (itsOK && track->GetTPCdEdxInfo(tpcdEdxInfo)){

      Bool_t isOK=(tpcdEdxInfo.GetNumberOfCrossedRows(0)>kMindEdxClustersRegion);
      isOK&=(tpcdEdxInfo.GetNumberOfCrossedRows(1)>kMindEdxClustersRegion);
      isOK&=(tpcdEdxInfo.GetNumberOfCrossedRows(2)>kMindEdxClustersRegion);
      isOK&=((tpcdEdxInfo.GetSignalMax(0)>0) && (tpcdEdxInfo.GetSignalMax(1)>0) && (tpcdEdxInfo.GetSignalMax(2)>0));
      isOK&=((tpcdEdxInfo.GetSignalTot(0)>0) && (tpcdEdxInfo.GetSignalTot(1)>0) && (tpcdEdxInfo.GetSignalTot(2)>0));
      isOK&=(itsOK||trdOK);      // stronger pile-up cut requiring ITS or TRD

      if (isOK) {
        (*fCacheTrackCounters)[6]+=1;         // Counter with accepted TPC dEdx info
        (*fCacheTrackdEdxRatio)[0]+=TMath::Log(tpcdEdxInfo.GetSignalMax(3));
        (*fCacheTrackdEdxRatio)[1]+=TMath::Log(tpcdEdxInfo.GetSignalTot(3));
        (*fCacheTrackdEdxRatio)[2]+=TMath::Log(tpcdEdxInfo.GetSignalMax(0)/tpcdEdxInfo.GetSignalTot(0));
        (*fCacheTrackdEdxRatio)[3]+=TMath::Log(tpcdEdxInfo.GetSignalMax(1)/tpcdEdxInfo.GetSignalTot(1));
        (*fCacheTrackdEdxRatio)[4]+=TMath::Log(tpcdEdxInfo.GetSignalMax(2)/tpcdEdxInfo.GetSignalTot(2));
        (*fCacheTrackdEdxRatio)[5]+=TMath::Log(tpcdEdxInfo.GetSignalMax(3)/tpcdEdxInfo.GetSignalTot(3));
        (*fCacheTrackdEdxRatio)[6]+=TMath::Log(tpcdEdxInfo.GetSignalMax(1)/tpcdEdxInfo.GetSignalMax(0));
        (*fCacheTrackdEdxRatio)[7]+=TMath::Log(tpcdEdxInfo.GetSignalMax(1)/tpcdEdxInfo.GetSignalMax(2));
        (*fCacheTrackdEdxRatio)[8]+=TMath::Log(tpcdEdxInfo.GetSignalTot(1)/tpcdEdxInfo.GetSignalTot(0));
        (*fCacheTrackdEdxRatio)[9]+=TMath::Log(tpcdEdxInfo.GetSignalTot(1)/tpcdEdxInfo.GetSignalTot(2));
      }
    }

    if (itsOK) {  // ITS
      (*fCacheTrackCounters)[0]++;
      (*fCacheTrackNcl)[0] += nclITS;
      (*fCacheTrackChi2)[0] += TMath::Min(TMath::Sqrt(chi2ITS),10.); // cutoff chi2 10
      (*fCacheTrackMatchEff)[2]+=trdOK;
      (*fCacheTrackMatchEff)[3]+=tofOK;
      (*fCacheTrackChi2)[4]+= (chi2TPC>0) ? TMath::Sqrt(chi2TPC):2; // TPC chi2 in case prolongation to ITS
      // long tracks properties
      if (nclITS>4){
        (*fCacheTrackCounters)[7]++;
        (*fCacheTrackNcl)[7] += nclITS;
        (*fCacheTrackChi2)[7]+=TMath::Min(TMath::Sqrt(chi2ITS),10.);
      }
    }
    if (trdOK) {// TRD    ///TODO - why chi2TRD could be smaller than 0?
      (*fCacheTrackCounters)[2]++;
      (*fCacheTrackNcl)[2] += nclTRD;
      (*fCacheTrackChi2)[2] += TMath::Sqrt(chi2TRD);
      (*fCacheTrackMatchEff)[0]+=itsOK;
      (*fCacheTrackChi2)[5]+= (chi2TPC>0) ? TMath::Sqrt(chi2TPC):2; // TPC chi2 in case prolongation to TRD
      if (nclTRD>80){
        (*fCacheTrackCounters)[8]++;
        (*fCacheTrackNcl)[8] += nclTRD;
        (*fCacheTrackChi2)[8]+=TMath::Min(TMath::Sqrt(chi2TRD),10.);
      }
    }
    if (tofOK) {  // TOF
      (*fCacheTrackCounters)[3]++;
      (*fCacheTrackNcl)[3] += 1;   // dummy for the moment
      (*fCacheTrackChi2)[3]+= 1;   //
    }
  } // end of track LOOP

  for (Int_t i=0; i<9; i++) if ((*fCacheTrackCounters)[i]>0) (*fCacheTrackNcl)[i]/=(*fCacheTrackCounters)[i];
  for (Int_t i=0; i<4; i++) if ((*fCacheTrackCounters)[i]>0) (*fCacheTrackChi2)[i]/=(*fCacheTrackCounters)[i];

  for (Int_t i=4; i<7; i++)  if ((*fCacheTrackCounters)[1]>0) (*fCacheTrackNcl)[i]/=(*fCacheTrackCounters)[1];
  if ((*fCacheTrackCounters)[6]>0){
    for (Int_t i=0; i<10; i++)   (*fCacheTrackdEdxRatio)[i]/=(*fCacheTrackCounters)[6];
  }
  //
  // conditional matching efficiency and chi2
  if ((*fCacheTrackCounters)[0]>0){
    (*fCacheTrackMatchEff)[2]/=(*fCacheTrackCounters)[0];  // TRD if ITS
    (*fCacheTrackMatchEff)[3]/=(*fCacheTrackCounters)[0];  // TOF if ITS
    (*fCacheTrackChi2)[4]/=(*fCacheTrackCounters)[0];
  }
  if ((*fCacheTrackCounters)[2]>0) {
    (*fCacheTrackMatchEff)[0]/=(*fCacheTrackCounters)[2];
    (*fCacheTrackChi2)[5]/=(*fCacheTrackCounters)[2];
  } //ITS if TRD
  (*fCacheTrackCounters)[9]=fEvent->GetNumberOfTracks();  // original number of ESDtracks
  return 1;
}

//
void AliESDtools::StreamEventVariables() {
  (*fStreamer) << "evVars" <<
    "gid=" << fGlobalID <<
    "trackCounters.=" << fCacheTrackCounters <<
    "trackTPCCountersZ.=" << fCacheTrackTPCCountersZ <<
    "trackdEdxRatio.=" << fCacheTrackdEdxRatio <<
    "trackNcl.=" << fCacheTrackNcl <<
    "trackChi2.=" << fCacheTrackChi2 <<
    "trackMatchEff.=" << fCacheTrackMatchEff <<
    "\n";
}

//
Int_t   AliESDtools::GetNearestTrack(Int_t indexTrk, Int_t paramType, AliTrackerBase *tracker){
  //
  const AliESDtrack *trackMatch = fEvent->GetTrack(indexTrk);
  //
  if ( !trackMatch->IsOn(0x10) || trackMatch->IsOn(0x1) ) {
    // check only TPC only tracks
    return -1;
  }
  //
  Int_t nTracks = fEvent->GetNumberOfTracks();
  const Double_t ktglCut = .05;
  const Double_t kqptCut = .99;
  const Double_t kAlphaCut = .5;
  //
  Double_t chi2Min = 10000;
  Int_t indexMin = -1;

  ULong_t trackStatusNearest = 0;
  ULong_t trackStatusProbe = trackMatch->GetStatus();

  AliESDtrack dbgTrk(*trackMatch);

  int status = 0;

  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    if (iTrack == indexTrk) {
      status = 2;
      continue;
    }
    const AliESDtrack *ptrack = fEvent->GetTrack(iTrack);
    if (ptrack==NULL) {
      status = -2;
      continue;
    }
    // status flags: 0x1 - ITSin, 0x10 - TPCin, 0x100 - TRDin, etc
    if (ptrack->IsOn(0x10)) {
      status = 4;
      continue;
    }
    if (ptrack->GetKinkIndex(0) < 0) continue;            // skip kink daughters
    const AliExternalTrackParam * track = nullptr;        //
    if (paramType == 0) track = ptrack;                   // Global track
    if (paramType == 1) track = ptrack->GetConstrainedParam();  // TPC only track at inner wall of TPC
    if (track == NULL) {
      status = -4;
      continue;
    }

    // first rough cuts
    // fP3 cut
    if ( TMath::Abs( track->GetTgl() - trackMatch->GetTgl() ) > ktglCut ) {
      status = 6;
      continue;
    }
    // fP4 cut
    if ( TMath::Abs( track->GetSigned1Pt() - trackMatch->GetSigned1Pt() ) > kqptCut) {
      status = 8;
      continue;
    }
    // fAlpha cut
    if ( TMath::Abs( track->GetAlpha() - trackMatch->GetAlpha() ) > kAlphaCut ) {
      status = 10;
      continue;
    }
    // calculate and extract track with smallest chi2 distance
    AliExternalTrackParam param(*track);
    if ( param.Rotate( trackMatch->GetAlpha() ) == kFALSE ) {
      status = 12;
      continue;
    }
    if (tracker) {
      if ( !tracker->PropagateTrackToBxByBz(&param, trackMatch->GetX(), trackMatch->GetMass(), 2.0, kFALSE, 0.8) ) {
        status = 14;
        continue;
      }
    }
    else {
      if ( !param.PropagateTo(trackMatch->GetX(), trackMatch->GetBz()) )
        continue;
    }
    Double_t chi2 = trackMatch->GetPredictedChi2(&param);
    // output debug information
    //
    if (chi2 < chi2Min){
      indexMin = iTrack;
      chi2Min = chi2;
      trackStatusNearest = ptrack->GetStatus();
    }
  }

  if (fStreamer) {
    (*fStreamer) << "debug" <<
      "status=" << status <<
      "\n";
  }

  if (fStreamer && indexMin > 0) {
    Double_t probeAlpha = trackMatch->GetAlpha();
    Double_t probeX     = trackMatch->GetX();
    AliESDtrack dbgTrkNearest(*(fEvent->GetTrack(indexMin)));
    (*fStreamer) << "nearestTrk" <<
      "gid=" << fGlobalID <<
      "probeAlpha=" << probeAlpha <<
      "probeX=" << probeX <<
      "indexTrack=" << indexTrk <<
      "indexNearest=" << indexMin <<
      "trackStatusProbe=" << trackStatusProbe <<
      "trackStatusNearest=" << trackStatusNearest <<
      "chi2=" << chi2Min <<
      "paramProbe.=" << &dbgTrk <<
      "paramNearest.=" << &dbgTrkNearest <<
      "\n";
  }

  return indexMin;

}


///
/// \param esdEvent   -
/// \param esdFriend  - in case ESD friend not avaliable - ITS tracks from vertex to be used
/// \param pcstream   - debug output
void AliESDtools::ProcessITSTPCmatchOut(AliESDEvent *const esdEvent, AliESDfriend *const esdFriend, TTreeStream *pcstream){
  //
  // Process ITS standalone tracks find match with closest TPC(or combined tracks) tracks
  // marian.ivanov@cern.ch
  // 0.) Init variables
  // 1.) GetTrack parameters at TPC inner wall
  // 2.) Match closest TPC  track  (STANDALONE/global) - chi2 match criteria
  //
  // Logic to be used in reco:
  // 1.) Find matching ITSalone->TPCalone
  // 2.) if (!TPCalone.FindClose(TPCother))  TPCalone.Addopt(ITSalone)
  // 3.) ff ((ITSalone.FindClose(Global)==0) CreateGlobaltrack
  const Double_t radiusMatch=84.;    // redius to propagate
  //
  const Double_t dFastPhiCut=0.2;        // 6 sigma (200 MeV) fast angular cut
  const Double_t dFastThetaCut=0.12;     // 6 sigma (200 MeV) fast angular cut
  const Double_t dFastPosPhiCut=0.06;    // 6 sigma (200 MeV) fast angular cut
  const Double_t dFastZCut=6;            // 6 sigma (200 MeV) fast  z difference cut
  const Double_t dFastPtCut=2.;          // 6 sigma (200 MeV) fast 1/pt cut
  const Double_t chi2Cut=100;            // chi2 matching cut
  //
  if (!esdFriend) return;  // not ITS standalone track
  if (esdFriend->TestSkipBit()) return; // friends tracks  not stored
  Int_t ntracks=esdEvent->GetNumberOfTracks();
  Float_t bz = esdEvent->GetMagneticField();
  //
  // 0.) Get parameters in reference radius TPC Inner wall
  //
  //
  TMatrixD vecPosR0(ntracks,6);   // possition and  momentum estimate at reference radius
  TMatrixD vecMomR0(ntracks,6);   //
  TMatrixD vecPosR1(ntracks,6);   // possition and  momentum estimate at reference radius TPC track
  TMatrixD vecMomR1(ntracks,6);   //
  Double_t xyz[3], pxyz[3];      //
  for (Int_t iTrack=0; iTrack<ntracks; iTrack++){
    AliESDtrack *track = esdEvent->GetTrack(iTrack);
    if(!track) continue;
    if (track->GetInnerParam()){
      const AliExternalTrackParam *trackTPC=track->GetInnerParam();
      trackTPC->GetXYZAt(radiusMatch,bz,xyz);
      trackTPC->GetPxPyPzAt(radiusMatch,bz,pxyz);
      for (Int_t i=0; i<3; i++){
        vecPosR1(iTrack,i)=xyz[i];
        vecMomR1(iTrack,i)=pxyz[i];
      }
      vecPosR1(iTrack,3)= TMath::ATan2(xyz[1],xyz[0]);    // phi pos angle
      vecMomR1(iTrack,3)= TMath::ATan2(pxyz[1],pxyz[0]);  // phi mom angle
      vecMomR1(iTrack,4)= trackTPC->GetSigned1Pt();;
      vecMomR1(iTrack,5)= trackTPC->GetTgl();;
    }
    AliESDfriendTrack* friendTrack=esdFriend->GetTrack(iTrack);
    if(!friendTrack) continue;
    if (friendTrack->GetITSOut()){
      const AliExternalTrackParam *trackITS=friendTrack->GetITSOut();
      trackITS->GetXYZAt(radiusMatch,bz,xyz);
      trackITS->GetPxPyPzAt(radiusMatch,bz,pxyz);
      for (Int_t i=0; i<3; i++){
        vecPosR0(iTrack,i)=xyz[i];
        vecMomR0(iTrack,i)=pxyz[i];
      }
      vecPosR0(iTrack,3)= TMath::ATan2(xyz[1],xyz[0]);
      vecMomR0(iTrack,3)= TMath::ATan2(pxyz[1],pxyz[0]);
      vecMomR0(iTrack,4)= trackITS->GetSigned1Pt();;
      vecMomR0(iTrack,5)= trackITS->GetTgl();;
    }
  }
  //
  // 1.) Find closest matching tracks, between the ITS standalone track
  // and  the all other tracks
  //  a.) caltegory  - All
  //  b.) category   - without ITS
  //
  //
  Int_t ntracksPropagated=0;
  AliExternalTrackParam extTrackDummy;
  AliESDtrack           esdTrackDummy;
  AliExternalTrackParam itsAtTPC;
  AliExternalTrackParam itsAtITSTPC;
  for (Int_t iTrack0=0; iTrack0<ntracks; iTrack0++){
    AliESDtrack *track0 = esdEvent->GetTrack(iTrack0);
    if(!track0) continue;
    if (track0->IsOn(AliVTrack::kTPCin)) continue;
    AliESDfriendTrack* friendTrack0=esdFriend->GetTrack(iTrack0);
    if (!friendTrack0) continue;
    //if (!track0->IsOn(AliVTrack::kITSpureSA)) continue;
    //if (!friendTrack0->GetITSOut()) continue;  // is there flag for ITS standalone?
    ntracksPropagated++;
    //
    // 2.) find clostest TPCtrack
    //     a.) all tracks
    Double_t minChi2All=10000000;
    Double_t minChi2TPC=10000000;
    Double_t minChi2TPCITS=10000000;
    Int_t indexAll=-1;
    Int_t indexTPC=-1;
    Int_t indexTPCITS=-1;
    Int_t ncandidates0=0; // n candidates - rough cut
    Int_t ncandidates1=0; // n candidates - rough + chi2 cut
    itsAtTPC=*(friendTrack0->GetITSOut());
    itsAtITSTPC=*(friendTrack0->GetITSOut());
    for (Int_t iTrack1=0; iTrack1<ntracks; iTrack1++){
      AliESDtrack *track1 = esdEvent->GetTrack(iTrack1);
      if(!track1) continue;
      if (!track1->IsOn(AliVTrack::kTPCin)) continue;
      // fast checks
      //
      if (TMath::Abs(vecPosR1(iTrack1,2)-vecPosR0(iTrack0,2))>dFastZCut) continue;
      if (TMath::Abs(vecPosR1(iTrack1,3)-vecPosR0(iTrack0,3))>dFastPosPhiCut) continue;
      if (TMath::Abs(vecMomR1(iTrack1,3)-vecMomR0(iTrack0,3))>dFastPhiCut) continue;
      if (TMath::Abs(vecMomR1(iTrack1,5)-vecMomR0(iTrack0,5))>dFastThetaCut) continue;
      if (TMath::Abs(vecMomR1(iTrack1,4)-vecMomR0(iTrack0,4))>dFastPtCut) continue;
      ncandidates0++;
      //
      const AliExternalTrackParam * param1= track1->GetInnerParam();
      if (!friendTrack0->GetITSOut()) continue;
      AliExternalTrackParam outerITS = *(friendTrack0->GetITSOut());
      if (!outerITS.Rotate(param1->GetAlpha())) continue;
      if (!outerITS.PropagateTo(param1->GetX(),bz)) continue; // assume track close to the TPC inner wall
      Double_t chi2 =  outerITS.GetPredictedChi2(param1);
      if (chi2>chi2Cut) continue;
      ncandidates1++;
      if (chi2<minChi2All){
        minChi2All=chi2;
        indexAll=iTrack1;
      }
      if (chi2<minChi2TPC && track1->IsOn(AliVTrack::kITSin)==0){
        minChi2TPC=chi2;
        indexTPC=iTrack1;
        itsAtTPC=outerITS;
      }
      if (chi2<minChi2TPCITS && track1->IsOn(AliVTrack::kITSin)){
        minChi2TPCITS=chi2;
        indexTPCITS=iTrack1;
        itsAtITSTPC=outerITS;
      }
    }
    //
    AliESDtrack * trackAll= (indexAll>=0)? esdEvent->GetTrack(indexAll):&esdTrackDummy;
    AliESDtrack * trackTPC= (indexTPC>=0)? esdEvent->GetTrack(indexTPC):&esdTrackDummy;
    AliESDtrack * trackTPCITS= (indexTPCITS>=0)? esdEvent->GetTrack(indexTPCITS):&esdTrackDummy;
    (*pcstream)<<"itsTPC"<<
                       "indexAll="<<indexAll<<          // index of closest track (chi2)
                       "indexTPC="<<indexTPC<<          // index of closest TPCalone tracks
                       "indexTPCITS="<<indexTPCITS<<    // index of closest cobined tracks
                       "ncandidates0="<<ncandidates0<<  // number of candidates
                       "ncandidates1="<<ncandidates1<<
                       //
                       "chi2All="<<minChi2All<<         // chi2 of closest  tracks
                       "chi2TPC="<<minChi2TPC<<
                       "chi2TPCITS="<<minChi2TPCITS<<
                       //
                       "track0.="<<track0<<             // ITS standalone tracks
                       "trackAll.="<<trackAll<<         // Closets other track
                       "trackTPC.="<<trackTPC<<         // Closest TPC only track
                       "trackTPCITS.="<<trackTPCITS<<   // closest combined track
                       //
                       "itsAtTPC.="<<&itsAtTPC<<        // ITS track parameters at the TPC alone track  frame
                       "itsAtITSTPC.="<<&itsAtITSTPC<<  // ITS track parameters at the TPC combeined track  frame
                       "\n";
  }
}



