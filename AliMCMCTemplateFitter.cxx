// Code by Martin Völkl:  martin.andreas.volkl@cern.ch

#include "AliMCMCTemplateFitter.h"
#include <TH1.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TMinuit.h>
#include <iostream>
#include <Math/SpecFuncMathMore.h>
#include <Math/QuantFuncMathMore.h>
#include <chrono>

ClassImp(AliMCMCTemplateFitter);

AliMCMCTemplateFitter::AliMCMCTemplateFitter(Int_t nMCFunctions, TH1D ** MCDistributions, TH1D * dataHistogram)
{
  // Constructor
  rd = new TRandom3(4); // set 0?
  TextOutput = kTRUE;
  fMakepjCorrelations = kTRUE;
  fMakeAjiCorrelations = kFALSE;
  fNIterations = 5000;
  fMaxCorrelation = 100;
  fBurnInRatio = 0.1;
  fSafetyFactor = 1.5;
  fFailedConvergences = 0;
  IsBurnIn = kFALSE;
  fNPar=nMCFunctions;
  fFitDiagsMC=MCDistributions;
  fFitDiagData=dataHistogram;
  
  fHasMonotonyPrior = new Bool_t[fNPar];
  fMonotonyPriorLowBin = new Int_t[fNPar];
  fMonotonyPriorHighBin = new Int_t[fNPar];
  for(int i=0;i<fNPar;i++) // Setup monotony priors
  {
      fHasMonotonyPrior[i] = kFALSE;
      fMonotonyPriorLowBin[i] = 0;
      fMonotonyPriorHighBin[i] = 0;
  }
  fNBins = dataHistogram->GetNbinsX();
  // Create state vector
  fAji = new double[fNPar*fNBins];
  fpj = new double[fNPar];
  fpjMarginals = new TH1D*[fNPar];
  //Set starting values for the parameters
  for(int j=0;j<fNPar;j++)
  {
      for(int i=0;i<fNBins;i++)
      {
          //fAji[j*fNBins+i]=TMath::Sqrt(MCDistributions[j]->GetBinContent(i+1))+0.1; // + 0.1 just for safety
          fAji[j*fNBins+i]=MCDistributions[j]->GetBinContent(i+1)+0.1; // + 0.1 just for safety
      }
      fpj[j]=fFitDiagData->Integral()/MCDistributions[j]->Integral()/double(fNPar); // Guess: Each template contributes equally
  }
  // Create Output Histograms
  fpjMarginals = new TH1D*[fNPar];
  for(int j=0;j<fNPar;j++)
  {
      fpjMarginals[j]=new TH1D(Form("p%dMarginal", j), "", 1000, 0., fFitDiagData->Integral()/MCDistributions[j]->Integral()*fSafetyFactor);
      fpjMarginals[j]->GetXaxis()->SetTitle(Form("p_{%d}",j));
      fpjMarginals[j]->GetYaxis()->SetTitle(Form("counts"));
  }
  fpj2DMarginals = new TH2D*[fNPar*fNPar];
  for(int j=0;j<fNPar;j++)
      for(int j2=0;j2<fNPar;j2++)
      {
          fpj2DMarginals[j*fNPar+j2] = new TH2D(Form("p%dp%dMarginal", j,j2), "", 200, 0., fFitDiagData->Integral()/MCDistributions[j]->Integral()*fSafetyFactor, 200, 0., fFitDiagData->Integral()/MCDistributions[j2]->Integral()*fSafetyFactor);
          fpj2DMarginals[j*fNPar+j2]->GetXaxis()->SetTitle(Form("p_{%d}",j));
          fpj2DMarginals[j*fNPar+j2]->GetYaxis()->SetTitle(Form("p_{%d}",j2));
      }
  fAjiMarginals = new TH2D*[fNPar];
  for(int j=0;j<fNPar;j++)
  {
      // This is t
      if(MCDistributions[j]->GetXaxis()->GetXbins()->GetSize() > 0)
          fAjiMarginals[j]=new TH2D(Form("A%diMarginal", j), "", MCDistributions[j]->GetNbinsX(), MCDistributions[j]->GetXaxis()->GetXbins()->GetArray(), 1000, 0., MCDistributions[j]->GetMaximum()+5. * TMath::Sqrt(MCDistributions[j]->GetMaximum()));
      else
          fAjiMarginals[j]=new TH2D(Form("A%diMarginal", j), "", MCDistributions[j]->GetNbinsX(), MCDistributions[j]->GetXaxis()->GetBinLowEdge(1), MCDistributions[j]->GetXaxis()->GetBinLowEdge(MCDistributions[j]->GetNbinsX()+1), 1000, 0., MCDistributions[j]->GetMaximum()+5. * TMath::Sqrt(MCDistributions[j]->GetMaximum()));
      fAjiMarginals[j]->GetXaxis()->SetTitle(Form("i axis"));
      fAjiMarginals[j]->GetYaxis()->SetTitle(Form("A_{ji}"));
  }
  
  
  // Initialize bookkeeping
  fAccrejTries=0;
  fAccrejAcc=0;
  fConsideredParameter=0;
  
  // Initialize minimizer
  minimizingObject = new TMinuit(1);
  minimizingObject->SetName("AliMCMCTemplateFitterPpjMaximizer");
  minimizingObject->SetPrintLevel(-1);
  SomewhereElseInAliMCMCTemplateFitter::fitter=this;  // Workaround
  minimizingObject->SetFCN(SomewhereElseInAliMCMCTemplateFitter::Outsourcedlogfpj);    
}

AliMCMCTemplateFitter::~AliMCMCTemplateFitter()
{
  // Destructor
  delete[] fAji;
  delete[] fpj;
  for(int j=0;j<fNPar;j++)
      delete fpjMarginals[j];
  delete[] fpjMarginals;
  for(int j=0;j<fNPar*fNPar;j++)
      delete fpj2DMarginals[j];
  delete[] fpj2DMarginals;
  delete minimizingObject;
  delete[] fHasMonotonyPrior;
  delete[] fMonotonyPriorLowBin;
  delete[] fMonotonyPriorHighBin;
}

void AliMCMCTemplateFitter::MakeLimitsFromPrior(Int_t j, Int_t i, Double_t &lower, Double_t &upper)
{
    // Important: The initial point must conform to the prior!
    lower = 0.;
    upper = 10000.;
    if(fHasMonotonyPrior[j])
    {
        if(i < fMonotonyPriorLowBin[j])
        {
            if(i==0)
                lower = 0.;
            else
                lower = fAji[j*fNBins+i-1];
            if(i<fNBins-1)
                upper=fAji[j*fNBins+i+1];
        }
        if(i > fMonotonyPriorHighBin[j])
        {
            if(i==fNBins-1)
                lower = 0.;
            else
                lower = fAji[j*fNBins+i+1];
            if(i>0)
                upper=fAji[j*fNBins+i-1];
        }
        if(i >= fMonotonyPriorLowBin[j] && i <= fMonotonyPriorHighBin[j]) // In peak region
        { // regular monotony unless next to peak or at peak
            Int_t PeakPosition=fMonotonyPriorLowBin[j];
            Double_t CurrentMaximum = fAji[j*fNBins+fMonotonyPriorLowBin[j]];
            for(int m=fMonotonyPriorLowBin[j];m<=fMonotonyPriorHighBin[j];m++)
                if(fAji[j*fNBins+m]>CurrentMaximum)
                {
                    CurrentMaximum = fAji[j*fNBins+m];
                    PeakPosition = m;
                }
            // Found the current peak at PeakPosition
            if(i<PeakPosition)
            {
                upper = fAji[j*fNBins+i+1];
                if(i == PeakPosition-1) upper = 10000.; // Can become new peak
                if(i==0) lower = 0.;
                else lower = fAji[j*fNBins+i-1];
            }
            if(i>PeakPosition)
            {
                upper = fAji[j*fNBins+i-1];
                if(i == PeakPosition+1) upper = 10000.; // Can become new peak
                if(i==fNBins-1) lower = 0.;
                else lower = fAji[j*fNBins+i+1];
            }
            if(i == PeakPosition)
            {
                upper = 10000.;
                if(i==0 || i==fNBins-1) lower = 0.;
                else lower = TMath::Min(fAji[j*fNBins+i+1], fAji[j*fNBins+i-1]);
            }
        }   
    }
    if(lower >= upper)
    {
        if(TMath::Abs(lower-upper)>1.e-10 && !IsBurnIn)
            cout << "monotony not given. lower= " << lower << " upper= " << upper << endl;
        double temp = lower;
        lower = upper;
        upper = temp;
    }
}

void AliMCMCTemplateFitter::Calulatelogfpj(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    // Calculate p(pj|others) to find maximum for acc-rej
    double x = par[0];  
    double fi=0.;
    double di=0.;
    f=0.;
    for(int i=0;i<fNBins;i++)
    {
        di = fFitDiagData->GetBinContent(i+1);
        fi=0.;
        for(int j=0;j<fNPar;j++)
            if(j != fConsideredParameter) fi+= fpj[j] * fAji[j*fNBins+i];
        fi+= x * fAji[fConsideredParameter*fNBins+i];
        if(fi > 1.e-10) // For now just ignore bins without content
            f += di * TMath::Log(fi) - fi;
    }
    f = -f; // uses Minimizer
}

double AliMCMCTemplateFitter::SamplePoissonLikelihoodInvGamma(Int_t n, Double_t a, Double_t b)
{
    double Xlow = TMath::Gamma(double(n+1), a) / 1.;
    double Xhigh = TMath::Gamma(double(n+1), b) / 1.;
    double x = Xlow + (Xhigh - Xlow) * rd->Rndm();
    double returnValue = ::ROOT::MathMore::gamma_quantile(x, n+1, 1.);
    Int_t iterations=0;
    if(isinf(returnValue))cout << "Value is infinite" << endl;
    while(isnan(returnValue) && iterations<10) // This is necessary, because the function does not always converge
    {
        fFailedConvergences++;
        iterations++;
        x = Xlow + (Xhigh - Xlow) * rd->Rndm();
        returnValue = ::ROOT::MathMore::gamma_quantile(x, double(n+1), 1.);
    }
    if(iterations==10) cout << "There was a problem with drawing the gamma quantile functions for n+1=" << n+1 << endl;
    return returnValue;
}


double AliMCMCTemplateFitter::SampleExpOverX(double a, double b)
{
    double integral1 = TMath::Log(TMath::Min(b,1.)) - TMath::Log(TMath::Min(a,1.));
    double integral2 = TMath::Exp(-TMath::Max(a,1.)) - TMath::Exp(-TMath::Max(b,1.));
    
    Bool_t acc= kFALSE;
    
    double localb,locala;
    double x=0.;
    
    while(!acc)
    {
        if(rd->Rndm() < integral1/(integral1+integral2)) // Sample below 1
        {
            locala = TMath::Min(a,1.);
            localb = TMath::Min(b,1.);
            x=TMath::Exp(rd->Rndm()*TMath::Log(localb/locala)+TMath::Log(locala));
            if(rd->Rndm()<TMath::Exp(-x)) // acceptance Ratio is exp(-x)
                acc=true;
        }
        else // Sample above 1
        {
            locala = TMath::Max(a,1.);
            localb = TMath::Max(b,1.);
            x=-TMath::Log(rd->Rndm()*(TMath::Exp(-localb)-TMath::Exp(-locala))+TMath::Exp(-locala));
            if(rd->Rndm()<1./x) // acceptance Ratio is exp(-x)
                acc=true;
        }
    }
    return x;
}

double AliMCMCTemplateFitter::SampleConditionalProbabilitypj(Int_t j)
{
    double pjmax = fFitDiagData->Integral()/fFitDiagsMC[j]->Integral()*fSafetyFactor;
    double pjymax = 0.;
    double di=0.;
    double logfrescaling=0.;
    double x,y,fy, fi;
    /*for(int i=0;i<fNBins;i++) // Quite broad estimate for the maximum - might lead to small acceptance rates
    {
        di = fFitDiagData->GetBinContent(i+1);
        if(IsBurnIn)
        {
            fi=0.;
            for(int jall=0;jall<fNPar;jall++)
                fi+= fpj[jall] * fAji[jall*fNBins+i];
            pjymax += di * TMath::Log(fi)-fi;  // relax condition for burn-in for more acceptance
        }
        else
            pjymax += di * (TMath::Log(di)-1.); // leads to problems if di = 0
    }*/
    
    Int_t temp=0.;
    fConsideredParameter=j;
    minimizingObject->mnparm(0, Form("pjparameter"), fpj[j], fpjMarginals[j]->GetXaxis()->GetXmax()/20., 0., fpjMarginals[j]->GetXaxis()->GetXmax(), temp);
    minimizingObject->Migrad();
    double maximMax = -minimizingObject->fAmin + 0.1; //
    Double_t MaxPoint=0.; Double_t MaxPointError=0.;
    minimizingObject->GetParameter(0, MaxPoint, MaxPointError);
    
    pjymax=maximMax;
    logfrescaling=pjymax;
    pjymax-=logfrescaling;
    pjymax = TMath::Exp(pjymax);
    //if(!IsBurnIn) cout << "pjymax = " << TMath::Exp(logfrescaling) << " compared to actual maximum of " << TMath::Exp(maximMax) << endl;
    
    // approximate width
    Int_t npar=1; Double_t * gin = 0x0; Double_t f; Double_t f2; Double_t * par= new Double_t[1];
    par[0] = MaxPoint;
    Calulatelogfpj(npar, gin, f, par, 0);
    double CompPoint = (MaxPoint* 99. +pjmax)/100.;
    if(TMath::Abs(CompPoint-MaxPoint) < 0.0001*pjmax)
        CompPoint = MaxPoint* 99./100.;
    par[0] = CompPoint;
    Calulatelogfpj(npar, gin, f2, par, 0);
    double alpha = (f2-f) /(CompPoint-MaxPoint)/(CompPoint-MaxPoint);
    //cout << "Width: " << 1./TMath::Sqrt(alpha) << endl;
    bool acc=false;
    while(!acc)
    {
        x=rd->Rndm()*pjmax;
        y=rd->Rndm()*pjymax;
        fy = 0.;
        
        for(int i=0;i<fNBins;i++)
        {
            fi = 0.;
            di = fFitDiagData->GetBinContent(i+1);
            fi += x*fAji[j*fNBins+i];
            for(int jother=0;jother<fNPar;jother++)
                if(jother != j) fi += fpj[jother] * fAji[jother*fNBins+i];
            if(fi > 1.e-10)
                fy += di * TMath::Log(fi) - fi;
        }
        fy-=logfrescaling;
        fy = TMath::Exp(fy);
        if(fy>pjymax && !IsBurnIn) cout << "Incorrect maximum in the accept-reject sampling of pj" << endl;
        if(y<fy)
            acc=true;
        fAccrejTries++;
    }
    fAccrejAcc++;
    return x;
}


double AliMCMCTemplateFitter::SampleConditionalPropabilityAjiajiIsMinusOne(Int_t j, Int_t i, Double_t a, Double_t b)
{  // Sample conditional B-B probability in range (a,b)
    Int_t aji = -1; // Changing aji is where a prior on the Aji could be applied
    if(a < 1.e-10) return 0; // Later do this with more cases
    if(b-a < 1.e-10) return a;
    Int_t di = fFitDiagData->GetBinContent(i+1);
    Double_t pj = fpj[j];
    Double_t fwoj = 0.;
    for(int jother=0;jother<fNPar;jother++)
        if(jother!=j) fwoj += fpj[jother]*fAji[jother*fNBins+i];
    Double_t aPrime = (1.+pj)*a; // cout << "aPrime: " << aPrime << " aji+1 = " << aji+1 << endl;
    Double_t bPrime = (1.+pj)*b;
    Double_t * Integral = new Double_t[di+1];
    Double_t Prefactor=0.;
    Double_t GammaA, GammaB;
    Double_t IntegralSum=0.;
    Double_t epsilon = 0.001; // The relative integral may be off by this amount
    Double_t IntegralMinusOne=0.;
    if(di==0) return SampleExpOverX(aPrime, bPrime)/(pj+1.);
    if(fwoj < 1.e-10) return SamplePoissonLikelihoodInvGamma(aji + di, aPrime, bPrime)/(pj+1.);
    IntegralMinusOne = TMath::Log((1.+pj))+ double(di)*TMath::Log(fwoj/pj);
    double minusoneintlow = 0.;
    double minusoneinthigh = 0.;
    if(aPrime < 50.) // hopefully, aPrime will always be small is aji=0
        minusoneintlow = -::ROOT::Math::expint(-aPrime);
    else
        minusoneintlow = 0.; // expint cannot handle large -x, but approaches 0 quickly
    if(bPrime < 50.) // hopefully, aPrime will always be small is aji=0
        minusoneinthigh = -::ROOT::Math::expint(-bPrime);
    else
        minusoneinthigh = 0.; // expint cannot handle large -x, but approaches 0 quickly
    if(minusoneintlow - minusoneinthigh > 1.e-100)
        IntegralMinusOne += TMath::Log(minusoneintlow - minusoneinthigh);
    else
        IntegralMinusOne = -1.e100;
    GammaA = 1. - TMath::Exp(-aPrime); //=TMath::Gamma(1., aPrime); // Important! TMath::Gamma is the LOWER incomplete Gamma function
    GammaB = 1. - TMath::Exp(-bPrime);//=TMath::Gamma(1., bPrime);
    if(GammaB-GammaA>0.)
        Integral[0] = TMath::Log(GammaB-GammaA) + TMath::Log(double(di)) + double(di-1)*TMath::Log(fwoj/pj); // Prefactor important to compare with l=0
    else
        Integral[0] = -1.e100;
    Double_t aSecondFactor = 0.;
    Double_t bSecondFactor = 0.;
    aSecondFactor = TMath::Exp(-aPrime);
    bSecondFactor = TMath::Exp(-bPrime);
    Prefactor = TMath::Log(double(di)) + double(di-1)*TMath::Log(fwoj/pj); // Normalization
    for(int l=2;l<=di;l++)
    {
        // Gamma function with x^(l-1) Gamma(x-2)
        aSecondFactor *= aPrime/double(l-1);
        bSecondFactor *= bPrime/double(l-1);
        Prefactor += TMath::Log(double(di-(l-1)) / double(l) / (pj+1.) / (fwoj/pj) * double(l-1)); // Use recursive formula
        GammaA = GammaA - aSecondFactor; // Also recursive formula for incomplete gamma function
        GammaB = GammaB - bSecondFactor;
        if(GammaB-GammaA>0) // Numerics occasionally gives negative integral it true difference is small
            Integral[l-1] = Prefactor + TMath::Log(GammaB-GammaA);
        else Integral[l-1]= -1.e200;
    }
    Int_t lmax = -1;
    Double_t maxInt = IntegralMinusOne;
    for(int l=0;l<=di-1;l++)
        if(Integral[l]>maxInt)
            maxInt=Integral[l];
    IntegralMinusOne -= maxInt;
    IntegralMinusOne = TMath::Exp(IntegralMinusOne);
    IntegralSum=IntegralMinusOne;
    for(int l=0;l<=di-1;l++)
    {
        Integral[l]-=maxInt;
        Integral[l] = TMath::Exp(Integral[l]);
        IntegralSum += Integral[l];
    }
    // Now find out from which of the terms one should draw the random number
    Double_t RDIntegral = IntegralSum * rd->Rndm();
    Double_t IntegralSumIteration = IntegralMinusOne;
    Int_t l = 0;
    while(IntegralSumIteration < RDIntegral && IntegralSumIteration<=IntegralSum) // Second is just for safety
        IntegralSumIteration += Integral[l++];
    l--;
    if(IntegralSumIteration>IntegralSum*(1.+epsilon))
        cout << "Unfortunately something has gone wrong when choosing the order of sampling function, the integral is off by " << (IntegralSumIteration/IntegralSum-1.)*100. << "% ( " << IntegralSumIteration << " vs " << IntegralSum << " )" << endl;
    // Now get random sample of the order aji+l between a and b (or aPrime,bPrime)
    Double_t RandomNumber = 0.;
    if(l>=0)
        RandomNumber = SamplePoissonLikelihoodInvGamma(l, aPrime, bPrime);
    else
        RandomNumber = SampleExpOverX(aPrime, bPrime);
    RandomNumber = RandomNumber/(pj+1.); // To transform back to Aji
    delete[] Integral;
    return RandomNumber;
}

double AliMCMCTemplateFitter::SampleConditionalProbabilityAji(Int_t j, Int_t i, Double_t a, Double_t b)
{  // Sample conditional B-B probability in range (a,b)
    Int_t aji = fFitDiagsMC[j]->GetBinContent(i+1) -1; // Changing aji, this is where a prior on the Aji could be applied
    if(aji<0) return SampleConditionalPropabilityAjiajiIsMinusOne(j, i, a, b); // Later do this with more cases
    if(b-a < 1.e-10) return a;
    Int_t di = fFitDiagData->GetBinContent(i+1);
    Double_t pj = fpj[j];
    Double_t fwoj = 0.;
    for(int jother=0;jother<fNPar;jother++)
        if(jother!=j) fwoj += fpj[jother]*fAji[jother*fNBins+i];
    Double_t aPrime = (1.+pj)*a; // cout << "aPrime: " << aPrime << " aji+1 = " << aji+1 << endl;
    Double_t bPrime = (1.+pj)*b;
    Double_t * Integral = new Double_t[di+1];
    Double_t Prefactor=0.;
    Double_t GammaA, GammaB;
    Double_t IntegralSum=0.;
    Double_t epsilon = 0.001; // The relative integral may be off by this amount
    if(fwoj<1.e-3) return SamplePoissonLikelihoodInvGamma(aji+di, aPrime, bPrime)/(pj+1.);
    // Do integrals for AjiPrime = Aji * (pj+1) , ratios of integrals stay the same, thus factor is ignored
    // Start with case l=0
    //Prefactor = TMath::Power(1./(pj+1.),aji+0) * TMath::Power(fwoj/pj,di);
    Prefactor = 1.; // Normalization
    // Now integral of x^aji exp(-x)
    //Double_t GammaNormalization = TMath::Gamma(double(aji+1)); // Gamma function is normalized. Need to undo this!
    GammaA = TMath::Gamma(double(aji+1), aPrime);//*GammaNormalization;  // Important! TMath::Gamma is the LOWER incomplete Gamma function
    GammaB = TMath::Gamma(double(aji+1), bPrime);//*GammaNormalization;
    if(GammaB-GammaA > 1.e-100)
        Integral[0] = TMath::Log(GammaB-GammaA);
    else
        Integral[0] = -1.e100;
    // Now we should have gamma(aji+1);
    Double_t aSecondFactor = 0.;
    Double_t bSecondFactor = 0.;
    if(aPrime>1.e-100){
    if(aji < 80)
        aSecondFactor = TMath::Exp(aji*TMath::Log(aPrime) -aPrime) / TMath::Gamma(aji+1);
    else // Stirling's approximation, crashes if a/bPrime=0
        aSecondFactor = TMath::Exp(double(aji)*TMath::Log(aPrime*TMath::E()/double(aji)) - aPrime - 0.5 *TMath::Log(2.*TMath::Pi() * double(aji)));
    }
    if(aji < 80)
    {
        bSecondFactor = TMath::Exp(aji*TMath::Log(bPrime) -bPrime) / TMath::Gamma(aji+1);
    }
    else
        bSecondFactor = TMath::Exp(double(aji)*TMath::Log(bPrime*TMath::E()/double(aji)) - bPrime - 0.5 *TMath::Log(2.*TMath::Pi() * double(aji)));
    //cout << "Test 2.5 aSecondFactor=" << aSecondFactor << " , bSecondFactor=" << bSecondFactor << endl;
    for(int l=1;l<=di;l++)
    {
        //ReducedGammaNormalization *= aji+l;
        aSecondFactor *= aPrime/double(aji+l);
        bSecondFactor *= bPrime/double(aji+l);
        Prefactor += TMath::Log(double(di-(l-1)) / double(l) / (pj+1.) / (fwoj/pj) * double(aji+l)); // Use recursive formula
        GammaA = GammaA - aSecondFactor; // Also recursive formula for incomplete gamma function
        GammaB = GammaB - bSecondFactor;
        if(GammaB-GammaA>0) // Numerics occasionally gives negative integral it true difference is small
            Integral[l] = Prefactor + TMath::Log(GammaB-GammaA);
        else Integral[l]= -1.e200;
    }
    Int_t lmax = 0;
    Double_t maxInt = -1.e100;
    //cout << "Integrals x " << di << endl;
    for(int l=0;l<=di;l++)
        if(Integral[l]>maxInt)
        {
            //cout << "Integral["<<l<<"] = " << Integral[l] << endl;
            maxInt=Integral[l];
            //cout << "end integral" << endl;
        }
    //cout << "Max Int= " << maxInt << endl;
    IntegralSum=0.;
    for(int l=0;l<=di;l++)
    {
        Integral[l]-=maxInt;
        Integral[l] = TMath::Exp(Integral[l]);
        IntegralSum += Integral[l];
    }
    // Now find out from which of the terms one should draw the random number
    Double_t RDIntegral = IntegralSum * rd->Rndm();
    Double_t IntegralSumIteration = 0.;
    Int_t l = 0;
    while(IntegralSumIteration < RDIntegral && IntegralSumIteration<=IntegralSum) // Second is just for safety
        IntegralSumIteration += Integral[l++];
    l--;
    //cout << "Integral is: " << l << " out of " << di <<endl;
    if(IntegralSumIteration>IntegralSum*(1.+epsilon))
        cout << "Unfortunately something has gone wrong when choosing the order of sampling function, the integral is off by " << (IntegralSumIteration/IntegralSum-1.)*100. << "% ( " << IntegralSumIteration << " vs " << IntegralSum << " )" << endl;
    
    // Now get random sample of the order aji+l between a and b (or aPrime,bPrime)
    Double_t RandomNumber = SamplePoissonLikelihoodInvGamma(aji+l, aPrime, bPrime);
    RandomNumber = RandomNumber/(pj+1.); // To transform back to Aji
    delete[] Integral;
    if(!(RandomNumber>a && RandomNumber<b)) cout << "Error: Random number " << RandomNumber << " is not between " << a << " and " << b << endl;
    return RandomNumber;
}

void AliMCMCTemplateFitter::DoOneIteration(void)
{
    double lowerLimit=0.; double upperLimit=0.;
    // Order for now: pj, then A1i, then A2i etc. Can be modified later
    for(int j=0;j<fNPar;j++)
    {
        fpj[j] = SampleConditionalProbabilitypj(j);
    }
    for(int j=0;j<fNPar;j++)
    {
        for(int i=0;i<fNBins;i++)
        {
            MakeLimitsFromPrior(j, i, lowerLimit, upperLimit);
            fAji[j*fNBins+i] = SampleConditionalProbabilityAji(j, i, lowerLimit, upperLimit); // later, limits can be set
        }
    }
}

TH1D * AliMCMCTemplateFitter::GetpjAutocorrelation(Int_t j, Int_t MaxLag)
{
    if(MaxLag>=fNIterations-1) MaxLag = fNIterations -2;
    double mean =0.;
    for(int i=0;i<fNIterations;i++) // Could have been filled with a different number of iterations
        mean+= fpjSteps[i*fNPar+j];
    mean /= double(fNIterations);
    double variance =0.;
    for(int i=0;i<fNIterations;i++)
        variance+= (fpjSteps[i*fNPar+j]-mean)*(fpjSteps[i*fNPar+j]-mean);
    variance /= double(fNIterations);
    
    TH1D * OutputCorr = new TH1D(Form("p%dCorrelationHist",j), "", MaxLag, 0.5, 0.5+double(MaxLag));
    OutputCorr->GetXaxis()->SetTitle("#DeltaIterations");
    OutputCorr->GetYaxis()->SetTitle("correlation coefficient");
    double corr=0.;
    for(int lag=1;lag<=MaxLag;lag++)
    {
        corr=0.;
        for(int i=0;i<fNIterations-lag;i++)
            corr += (fpjSteps[i*fNPar+j]-mean)*(fpjSteps[(i+lag)*fNPar+j]-mean)/variance;
        corr /= double(fNIterations-lag);
        OutputCorr->SetBinContent(lag,corr);
    }
    return OutputCorr;
}


TH1D * AliMCMCTemplateFitter::GetAjiAutocorrelation(Int_t j, Int_t i, Int_t MaxLag)
{
    if(MaxLag>=fNIterations-1) MaxLag = fNIterations -2;
    double mean =0.;
    for(int k=0;k<fNIterations;k++) // Could have been filled with a different number of iterations
        mean+= fAjiSteps[k*(fNPar*fNBins)+i*fNPar+j];
    mean /= double(fNIterations);
    double variance =0.;
    for(int k=0;k<fNIterations;k++)
        variance+= (fAjiSteps[k*(fNPar*fNBins)+i*fNPar+j]-mean)*(fAjiSteps[k*(fNPar*fNBins)+i*fNPar+j]-mean);
    variance /= double(fNIterations);
    TH1D * OutputCorr = new TH1D(Form("A%d%dCorrelationHist",j,i), "", MaxLag, 0.5, 0.5+double(MaxLag));
    OutputCorr->GetXaxis()->SetTitle("#DeltaIterations");
    OutputCorr->GetYaxis()->SetTitle("correlation coefficient");
    double corr=0.;
    // definition: If the variance is zero, the correlation is just 0
    for(int lag=1;lag<=MaxLag;lag++)
    {
        corr=0.;
        if(variance>0.001)
        for(int k=0;k<fNIterations-lag;k++)
            corr += (fAjiSteps[k*(fNPar*fNBins)+i*fNPar+j]-mean)*(fAjiSteps[(k+lag)*(fNPar*fNBins)+i*fNPar+j]-mean)/variance;
        corr /= double(fNIterations-lag);
        OutputCorr->SetBinContent(lag,corr);
    }
    return OutputCorr;
}

void AliMCMCTemplateFitter::PrepareCorrelations(void)
{
    fpjSteps = new Double_t[fNIterations*fNPar];
    if(fMakeAjiCorrelations) fAjiSteps = new Double_t[fNIterations*fNPar*fNBins];
}

void AliMCMCTemplateFitter::FillStatisticsObjects(Int_t iter)
{
    for(int j=0;j<fNPar;j++)
            fpjMarginals[j]->Fill(fpj[j]);
    for(int j=0;j<fNPar;j++)
      for(int j2=0;j2<fNPar;j2++)
          fpj2DMarginals[j*fNPar+j2]->Fill(fpj[j],fpj[j2]);
    for(int j=0;j<fNPar;j++)
        for(int i=0;i<fNBins;i++)
            fAjiMarginals[j]->Fill(fAjiMarginals[j]->GetXaxis()->GetBinCenter(i+1), fAji[j*fNBins+i]);
    
      if(fMakepjCorrelations)
        for(int j=0;j<fNPar;j++)
            fpjSteps[iter*fNPar+j] = fpj[j];
    if(fMakeAjiCorrelations)
        for(int j=0;j<fNPar;j++)
            for(int k=0;k<fNBins;k++)
            fAjiSteps[iter*(fNPar*fNBins)+k*fNPar+j] = fAji[j*fNBins+k];
}

void AliMCMCTemplateFitter::Fit(Int_t Iterations)
{
    double timediff=0.;
    auto t1 = std::chrono::high_resolution_clock::now();
    auto t2 = std::chrono::high_resolution_clock::now();
    fNIterations = Iterations;
    PrepareCorrelations();
    Int_t BurnInIterations = Int_t(fNIterations * fBurnInRatio);
    if(TextOutput) cout << "Starting burn-in with " << BurnInIterations << " iterations." << endl;
    IsBurnIn = kTRUE;
    t1 = std::chrono::high_resolution_clock::now();
    for(int iter=0;iter<BurnInIterations;iter++)
    {
        if(TextOutput) if((iter*100)/BurnInIterations>((iter-1)*100)/BurnInIterations) cout << "\r" << (iter*100)/BurnInIterations << "\%  " << flush;
        DoOneIteration();
    }
    t2 = std::chrono::high_resolution_clock::now();
    timediff = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    if(TextOutput) cout << "\rBurn-in done after " << int(timediff/1.e6) << "s. Average time per iteration: " << timediff/double(BurnInIterations) << " microseconds." << endl;
    if(TextOutput) cout << "Starting main fit with " << fNIterations<< " iterations." << endl;
    IsBurnIn = kFALSE;
    t1 = std::chrono::high_resolution_clock::now();
    fAccrejAcc=0;
    fAccrejTries=0;
    fFailedConvergences = 0; // Reset, only failed convergences in main fit are relevant
    for(int iter=0;iter<fNIterations;iter++)
    {
        if(TextOutput) if((iter*100)/fNIterations>((iter-1)*100)/fNIterations) cout << "\r" << (iter*100)/fNIterations << "\%  " << flush;
        DoOneIteration();
        FillStatisticsObjects(iter);
    }
    t2 = std::chrono::high_resolution_clock::now();
    timediff = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    if(TextOutput) cout << "\rMain fit done after " << int(timediff/1.e6) << "s. Average time per iteration: " << timediff/double(fNIterations) << " microseconds." << endl;
    if(TextOutput) cout << "Failed gamma_quantile calculation ratio: " << double(fFailedConvergences)/double(fNIterations) * 100. << "%" << endl;
    if(TextOutput) cout << "Acceptance ratio for Acc-Rej steps: " << double(fAccrejAcc)/double(fAccrejTries) * 100. << "%" << endl;
    
}

TH2D * AliMCMCTemplateFitter::GetCredibleArea(Int_t j, Int_t i, double probability, Bool_t smooth) // To be drawn as contour
{
    if(probability>1.) probability = 0.95;
    TH2D * contour = (TH2D*) fpj2DMarginals[j*fNPar+i]->Clone(Form("%s-contour", fpj2DMarginals[j*fNPar+i]->GetName()));
    int nx = contour->GetNbinsX();
    int ny = contour->GetNbinsY();
    double * histentries = new double[nx*ny];
    if(smooth) contour->Smooth();
    for(int i=0;i<nx;i++)
        for(int j=0;j<ny;j++)
            histentries[i*ny+j] = contour->GetBinContent(i+1,j+1);
    int * indices = new int[nx*ny];
    TMath::Sort(nx*ny,histentries, indices);
    Double_t sum = 0.;
    for(int i=0;i<nx*ny;i++)
        sum += histentries[indices[i]];
    Double_t sumIter = 0.;
    Int_t ele = 0;
    while(sumIter < sum * probability)
        sumIter += histentries[indices[ele++]];
    double * levels = new double[1];
    levels[0] = (histentries[indices[ele-1]] + histentries[indices[ele]]) / 2.;
    contour->SetContour(1, levels);
    contour->SetLineColor(kBlack);
    contour->SetLineWidth(2);
    return contour;
}
