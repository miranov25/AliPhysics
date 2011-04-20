#ifndef AliOADBTriggerAnalysis_H
#define AliOADBTriggerAnalysis_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     OADB container for filling scheme information (BX ids, name ...)
//     Author: Michele Floris, CERN
//-------------------------------------------------------------------------

#include <TNamed.h>
#include "TMap.h"
#include "TObjString.h"


class AliOADBTriggerAnalysis : public TNamed {

 public :
  AliOADBTriggerAnalysis();
  AliOADBTriggerAnalysis(char* name);
  virtual ~AliOADBTriggerAnalysis();
  //  void Init();
  
  // Getters
  Float_t GetZDCCutRefSumCorr()     { return fZDCCutRefSumCorr;     }      
  Float_t GetZDCCutRefDeltaCorr()   { return fZDCCutRefDeltaCorr;   }   
  Float_t GetZDCCutSigmaSumCorr()   { return fZDCCutSigmaSumCorr;   }   
  Float_t GetZDCCutSigmaDeltaCorr() { return fZDCCutSigmaDeltaCorr; }  
  // Setters
  void SetZDCCorrParameters(Float_t sumCorr, Float_t deltaCorr, Float_t sigmaSumCorr, Float_t sigmaDeltaCorr) 
  { fZDCCutRefSumCorr = sumCorr; fZDCCutRefDeltaCorr = deltaCorr; fZDCCutSigmaSumCorr = sigmaSumCorr; fZDCCutSigmaDeltaCorr = sigmaDeltaCorr;}
  // Browse
  virtual Bool_t	IsFolder() const { return kTRUE; }
  void Browse(TBrowser *b);
  // Print
  virtual void	Print(Option_t* option = "") const;

 private :

  Float_t fZDCCutRefSumCorr;      // Corrected ZDC time cut configuration
  Float_t fZDCCutRefDeltaCorr;    // Corrected ZDC time cut configuration
  Float_t fZDCCutSigmaSumCorr;    // Corrected ZDC time cut configuration
  Float_t fZDCCutSigmaDeltaCorr;  // Corrected ZDC time cut configuration  

  AliOADBTriggerAnalysis(const AliOADBTriggerAnalysis& cont);  // not implemented
  AliOADBTriggerAnalysis& operator=(const AliOADBTriggerAnalysis& cont); // not implemented


  ClassDef(AliOADBTriggerAnalysis, 1);
};

#endif
