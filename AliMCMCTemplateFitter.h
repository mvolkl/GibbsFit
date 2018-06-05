#ifndef ALIMCMCTEMPLATEFITTER_H
#define ALIMCMCTEMPLATEFITTER_H


// Code by Martin Völkl:  martin.andreas.volkl@cern.ch
#include "TObject.h"

class TH1;
class TH1D;
class TH2D;
class TRandom3;
class AliMCMCTemplateFitter;
class TMinuit;


class AliMCMCTemplateFitter
{
public:
  AliMCMCTemplateFitter(Int_t nMCFunctions, TH1D ** MCDistributions, TH1D * dataHistogram); 
  ~AliMCMCTemplateFitter();
  void DoOneIteration(void);
  void Fit(Int_t Iterations);
  TH1D * ReturnpjMarginal(Int_t j){return fpjMarginals[j];}
  TH2D * ReturnAjiMarginals(Int_t j){return fAjiMarginals[j];}
  TH2D * ReturnpjMarginal(Int_t j, Int_t j2){return fpj2DMarginals[j*fNPar+j2];}
  void Calulatelogfpj(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  //void EstimateIfConverged(); // compare to autocorrelation
  void SetPreparepjAutocorrelation(Bool_t PrepareAutocorrelation = kTRUE){fMakepjCorrelations = PrepareAutocorrelation;}
  TH1D * GetpjAutocorrelation(Int_t j, Int_t MaxLag);
  TH1D * GetAjiAutocorrelation(Int_t j, Int_t i, Int_t MaxLag);
  TH2D * GetCredibleArea(Int_t j, Int_t i, double probability = 0.95, Bool_t smooth = kTRUE); // To be drawn as contour
  void SetPrepareAjiAutocorrelation(Bool_t PrepareAutocorrelation = kTRUE){fMakeAjiCorrelations = PrepareAutocorrelation;}
  void SetMonotonyPrior(Int_t j, Int_t binBelowPeak, Int_t binAbovePeak){fHasMonotonyPrior[j]=kTRUE; fMonotonyPriorLowBin[j] = binBelowPeak; fMonotonyPriorHighBin[j]=binAbovePeak;}
  void SetBurnInRatio(Double_t BurnInRatio = 0.1){fBurnInRatio = BurnInRatio;}
  
private:
  double SamplePoissonLikelihoodInvGamma(Int_t n, Double_t a, Double_t b);
  double SampleConditionalProbabilityAji(Int_t j, Int_t i, Double_t a, Double_t b);
  double SampleConditionalPropabilityAjiajiIsMinusOne(Int_t j, Int_t i, Double_t a, Double_t b);
  double SampleExpOverX(double a, double b);
  double SampleConditionalProbabilitypj(Int_t j);
  void MakeLimitsFromPrior(Int_t j, Int_t i, Double_t &lower, Double_t &upper);
  void FillStatisticsObjects(Int_t iter);
  void PrepareCorrelations(void);
  Int_t fNPar; // Number of free parameters (number of templates for fitting)
  TH1D ** fFitDiagsMC; // The MC templates
  TH1D * fFitDiagData; // Data input histogram
  double * fAji; // consecutive arrays of the Aji with fixed j
  double * fpj;
  Int_t fNBins;
  Int_t fNIterations;
  Double_t fBurnInRatio;
  Bool_t * fHasMonotonyPrior;
  Int_t * fMonotonyPriorLowBin;
  Int_t * fMonotonyPriorHighBin;
  TRandom3 * rd;
  Double_t fSafetyFactor; // How much can pj go beyond integral(pj*Aji)=integral(data)
  TH1D ** fpjMarginals; // 1d Marginals of the pj for output
  TH2D ** fAjiMarginals; // 1d Aji magrinals for each source for output
  TH2D ** fpj2DMarginals; // 2d Marginals of the pj for output
  Bool_t TextOutput;
  Bool_t IsBurnIn;
  Bool_t fMakepjCorrelations;
  Bool_t fMakeAjiCorrelations;
  Int_t fAccrejTries;
  Int_t fAccrejAcc;
  Int_t fConsideredParameter;
  Int_t fMaxCorrelation;
  Int_t fFailedConvergences;
  Double_t * fpjSteps; // For the correlations
  Double_t * fAjiSteps; // For the correlations
  TMinuit * minimizingObject;
};


namespace SomewhereElseInAliMCMCTemplateFitter
{
  AliMCMCTemplateFitter * fitter;
  void Outsourcedlogfpj(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
  {
    fitter->Calulatelogfpj(npar, gin, f, par, iflag);
  }
}

#endif
