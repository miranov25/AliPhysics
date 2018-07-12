#ifndef ALIESDTOOLS_H
#define ALIESDTOOLS_H

class AliTrackerBase;

class AliESDtools : public TNamed
{
  public:
  AliESDtools();
  AliESDtools(const AliESDtools &tools);
  AliESDtools& operator=(const AliESDtools &rhs);
  void Init(TTree* tree);
  void DeInit();
  /// caching
  Int_t  CacheTPCEventInformation();
  Int_t CalculateEventVariables();
  void StreamEventVariables();
  void TPCVertexFit(TH1F *hisVertex);
  Int_t  GetNearestTrack(Int_t indexTrk, Int_t probeType, Int_t paramType, AliTrackerBase *tracker = 0x0);
  void   ProcessITSTPCmatchOut(AliESDEvent *const esdEvent, AliESDfriend *const esdFriend, TTreeStream *pcstream);
  // static functions for querying in TTree formula
  static Int_t    SCalculateEventVariables(Int_t entry){fgInstance->fESDtree->GetEntry(entry); return fgInstance->CalculateEventVariables();}
  static Double_t GetTrackCounters(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackCounters)[index];}
  static Double_t GetTrackTPCCountersZ(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackTPCCountersZ)[index];}
  static Double_t GetTrackdEdxRatio(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackdEdxRatio)[index];}
  static Double_t GetTrackNcl(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackNcl)[index];}
  static Double_t GetTrackChi2(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackChi2)[index];}
  static Double_t GetTrackMatchEff(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackMatchEff)[index];}
  static Double_t GetMeanHisTPCVertexA(){return fgInstance->fHisTPCVertexA->GetMean();}
  static Double_t GetMeanHisTPCVertexC(){return fgInstance->fHisTPCVertexC->GetMean();}
  //
  Int_t fVerbose;
  TTree *fESDtree;
  AliESDEvent * fEvent;
  TH1F *fHisTPCVertexA;
  TH1F *fHisTPCVertexC;
  TH1F *fHisTPCVertex;
  TH1F *fHisTPCVertexACut;
  TH1F *fHisTPCVertexCCut;
  TH1F             * fHistPhiTPCcounterA;         // helper histogram phi counteres
  TH1F             * fHistPhiTPCcounterC;         // helper histogram phi counters
  TH1F             * fHistPhiTPCcounterAITS;      // helper histogram phi counters
  TH1F             * fHistPhiTPCcounterCITS;      // helper histogram phi counters
  TH1F             * fHistPhiITScounterA;         // helper histogram phi counters
  TH1F             * fHistPhiITScounterC;         // helper histogram phi counters
  TVectorF         * fCacheTrackCounters;         // track counter
  TVectorF         * fCacheTrackTPCCountersZ;     // track counter with DCA z cut
  TVectorF         * fCacheTrackdEdxRatio;        // dedx info counter
  TVectorF         * fCacheTrackNcl;              // ncl counter
  TVectorF         * fCacheTrackChi2;             // chi2 counter
  TVectorF         * fCacheTrackMatchEff;         // matchEff counter
  TGraph           * fLumiGraph;                  // graph for the interaction rate info for a run
  ULong64_t fGlobalID;                            // global event ID
  //
  TTreeSRedirector * fStreamer;                   // streamer
  //
  static AliESDtools* fgInstance;                /// instance of the tool
  ClassDef(AliESDtools, 1) 
};

#endif
