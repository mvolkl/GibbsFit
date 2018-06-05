#include "MartinsStyle.h"
#include "AliMCMCTemplateFitter.h"
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

void SetHistogram(TH1D * in, int color)
{
    in->SetLineWidth(2);
    in->SetLineColor(color);
    in->SetMarkerColor(color);
}

double SamplePoissonLikelihoodInvGamma(Int_t n, Double_t a, Double_t b, TRandom3 * rd)
{
    double norm = TMath::Gamma(double(n+1));
    double Xlow = TMath::Gamma(double(n+1), a) / 1.;
    double Xhigh = TMath::Gamma(double(n+1), b) / 1.;
    double x = Xlow + (Xhigh - Xlow) * rd->Rndm();
    return ::ROOT::Math::gamma_quantile(x, double(n+1), 1.);
}

double SampleConditionalProbabilityAji(Int_t aji, Int_t di, Double_t fwoj, Double_t pj, Double_t a, Double_t b, TRandom3 * rd)
{  // Sample conditional B-B probability in range (a,b)
    Double_t aPrime = (1.+pj)*a;
    Double_t bPrime = (1.+pj)*b;
    Double_t * Integral = new Double_t[di+1];
    Double_t Prefactor=0.;
    Double_t GammaA, GammaB;
    Double_t IntegralSum=0.;
    Double_t epsilon = 0.001; // The relative integral may be off by this amount
    
    // Do integrals for AjiPrime = Aji * (pj+1) , ratios of integrals stay the same, thus factor is ignored
    // Start with case l=0
    //Prefactor = TMath::Power(1./(pj+1.),aji+0) * TMath::Power(fwoj/pj,di);
    Prefactor = 1.; // Normalization
    // Now integral of x^aji exp(-x)
    Double_t GammaNormalization = TMath::Gamma(double(aji+1)); // Gamma function is normalized. Need to undo this!
    GammaA = TMath::Gamma(double(aji+1), aPrime)*GammaNormalization;  // Important! TMath::Gamma is the LOWER incomplete Gamma function
    GammaB = TMath::Gamma(double(aji+1), bPrime)*GammaNormalization;
    Integral[0] = Prefactor * (GammaB-GammaA);
    IntegralSum=Integral[0];
    for(int l=1;l<=di;l++)
    {
        //Prefactor = TMath::Binomial(di,l) * TMath::Power(1./(pj+1.),aji+l) * TMath::Power(fwoj/pj,di-l);
        Prefactor = Prefactor * double(di-(l-1)) / double(l) / (pj+1.) / (fwoj/pj); // Use recursive formula
        GammaA = double(aji+l)*GammaA - TMath::Power(aPrime, aji+l) * TMath::Exp(-aPrime); // Also recursive formula for incomplete gamma function
        GammaB = double(aji+l)*GammaB - TMath::Power(bPrime, aji+l) * TMath::Exp(-bPrime);
        Integral[l] = Prefactor * (GammaB-GammaA);
        IntegralSum += Integral[l];
        //cout << "Integral["<<l<<"] : " << Integral[l] << " (GammaB-GammaA) = (" << GammaB << " - " << GammaA << " = " <<(GammaB-GammaA) << endl;
    }
    // Now find out from which of the terms one should draw the random number
    Double_t RDIntegral = IntegralSum * rd->Rndm();
    Double_t IntegralSumIteration = 0.;
    Int_t l = 0;
    while(IntegralSumIteration < RDIntegral && IntegralSumIteration<=IntegralSum) // Second is just for safety
        IntegralSumIteration += Integral[l++];
    l--;
    if(IntegralSumIteration>IntegralSum*(1.+epsilon))
        cout << "Unfortunately something has gone wrong when choosing the order of the sampling function, the integral is off by " << (IntegralSumIteration/IntegralSum-1.)*100. << "% ( " << IntegralSumIteration << " vs " << IntegralSum << " )" << endl;
    
    // Now get random sample of the order aji+l between a and b (or aPrime,bPrime)
    Double_t RandomNumber = SamplePoissonLikelihoodInvGamma(aji+l, aPrime, bPrime, rd);
    RandomNumber = RandomNumber/(pj+1.); // To transform back to Aji
    delete[] Integral;
    return RandomNumber;
}

int p1acc;
int p1all;

void DoGibbsStep2Sources(TH1D * t1, TH1D * t2, TH1D * data, Double_t * p, TH1D * A1i, TH1D * A2i, TRandom3 * rd)
{
    double safetyFactor = 5./TMath::Sqrt(data->Integral()); // 3 sigma
    double p1max = data->Integral()/t1->Integral()*(1.+safetyFactor);
    double p2max = data->Integral()/t2->Integral()*(1.+safetyFactor);
    double p1ymax = 0.;
    double p2ymax = 0.;
    double logfrescaling=0.;
    double di=0.;
    for(int i=1;i<=3;i++)
    {
        di = data->GetBinContent(i);
        p1ymax += di * (TMath::Log(di)-1.);
        p2ymax += di * (TMath::Log(di)-1.);
    }
    logfrescaling=p1ymax; // Log gets too large
    p1ymax -= logfrescaling;
    p2ymax -= logfrescaling;
    p1ymax = TMath::Exp(p1ymax);
    p2ymax = TMath::Exp(p2ymax);
    // Now sample the pj with accept-reject
    bool acc=false;
    double x,y,fy;
    while(!acc)
    {
        x=rd->Rndm()*p1max;
        y=rd->Rndm()*p1ymax;
        fy = 0.;
        for(int i=1;i<=3;i++)
        {
            di = data->GetBinContent(i);
            fy += di * TMath::Log(x*A1i->GetBinContent(i)+p[1]*A2i->GetBinContent(i)) - (x*A1i->GetBinContent(i)+p[1]*A2i->GetBinContent(i));
        }
        fy-=logfrescaling;
        fy = TMath::Exp(fy);
        if(fy>p1ymax) cout << "Incorrect maximum in the accept-reject sampling of p1" << endl;
        if(y<fy)
            acc=true;
        p1all++;
        //cout << "choose p1 between 0 and " << p1max << ": " << x << endl;
    }
    p1acc++;
    p[0]=x;
    acc=false;
    while(!acc)
    {
        x=rd->Rndm()*p2max;
        y=rd->Rndm()*p2ymax;
        fy = 0.;
        for(int i=1;i<=3;i++)
        {
            di = data->GetBinContent(i);
            fy += di * TMath::Log(p[0]*A1i->GetBinContent(i)+x*A2i->GetBinContent(i)) - (p[0]*A1i->GetBinContent(i)+x*A2i->GetBinContent(i));
        }
        fy-=logfrescaling;
        fy = TMath::Exp(fy);
        if(fy>p2ymax) cout << "Incorrect maximum in the accept-reject sampling of p2" << endl;
        if(y<fy)
            acc=true;
    }
    p[1]=x;
    double fwoj=0.;
    for(int i=1;i<=3;i++) // sample A1i
    {
        di = data->GetBinContent(i);
        fwoj = p[1]*A2i->GetBinContent(i);
         //cout << "calling SampleConditionalPropabilityAji(" <<  t1->GetBinContent(i) << ", " << di << ", " << fwoj <<", "<< p[0] << ", 0., 10., rd)" << endl;
        A1i->SetBinContent(i, SampleConditionalProbabilityAji(t1->GetBinContent(i), di, fwoj, p[0], 0., 1000., rd));
    }
    for(int i=1;i<=3;i++) // sample A2i
    {
        di = data->GetBinContent(i);
        fwoj = p[0]*A1i->GetBinContent(i);
        A2i->SetBinContent(i, SampleConditionalProbabilityAji(t2->GetBinContent(i), di, fwoj, p[1], 0., 1000., rd));
    }
}   

TH1D * A1iProposal;
TH1D * A2iProposal;
double * pprop;

Bool_t DoMetropolisStep2Sources(TH1D * t1, TH1D * t2, TH1D * data, Double_t * p, TH1D * A1i, TH1D * A2i, TRandom3 * rd) // Returns true is step was accepted
{
    double StepSizeFactor = 0.5;
    double StepSizep1 = 1. * StepSizeFactor;
    double StepSizep2 = 1. * StepSizeFactor;
    double CurrentLogProb = 0.;
    for(int i=1;i<=3;i++)
    {
        CurrentLogProb += data->GetBinContent(i)*TMath::Log(p[0]*A1i->GetBinContent(i)+p[1]*A2i->GetBinContent(i)) -(p[0]*A1i->GetBinContent(i)+p[1]*A2i->GetBinContent(i))
            + t1->GetBinContent(i)*TMath::Log(A1i->GetBinContent(i)) -(A1i->GetBinContent(i)) +t2->GetBinContent(i)*TMath::Log(A2i->GetBinContent(i)) -(A2i->GetBinContent(i));
    }
    // Make Proposal
    pprop[0]=  rd->Gaus(p[0], StepSizep1);
    pprop[1]=  rd->Gaus(p[1], StepSizep2);
    
    for(int i=1;i<=3;i++)
    {
        A1iProposal->SetBinContent(i, rd->Gaus(A1i->GetBinContent(i), StepSizeFactor*TMath::Sqrt(t1->GetBinContent(i))));
        A2iProposal->SetBinContent(i, rd->Gaus(A2i->GetBinContent(i), StepSizeFactor*TMath::Sqrt(t2->GetBinContent(i))));
    }
    
    double ProposalLogProb = 0.;
    // Check if any parameters are negative
    Bool_t finiteProbability = kTRUE;
    for(int i=1;i<=3;i++)
        if(A1iProposal->GetBinContent(i) < 0. || A2iProposal->GetBinContent(i)<0.) finiteProbability=kFALSE;
    if(pprop[0]<0. || pprop[1]<0.) finiteProbability=kFALSE;
    if(finiteProbability)
    for(int i=1;i<=3;i++)
    {
        ProposalLogProb += data->GetBinContent(i)*TMath::Log(pprop[0]*A1iProposal->GetBinContent(i)+pprop[1]*A2iProposal->GetBinContent(i)) -(pprop[0]*A1iProposal->GetBinContent(i)+pprop[1]*A2iProposal->GetBinContent(i))
            + t1->GetBinContent(i)*TMath::Log(A1iProposal->GetBinContent(i)) -(A1iProposal->GetBinContent(i)) +t2->GetBinContent(i)*TMath::Log(A2iProposal->GetBinContent(i)) -(A2iProposal->GetBinContent(i));
    }
    double AcceptanceRatio = 0.;
    if(finiteProbability)
        AcceptanceRatio = TMath::Exp(ProposalLogProb-CurrentLogProb);
    else
        AcceptanceRatio=0.;
    
    if(rd->Rndm()<AcceptanceRatio)
    {
        p[0]=pprop[0];
        p[1]=pprop[1];
        for(int i=1;i<=3;i++)
        {
            A1i->SetBinContent(i, A1iProposal->GetBinContent(i));
            A2i->SetBinContent(i, A2iProposal->GetBinContent(i));
        }
        return kTRUE;
    }
    return kFALSE;
}


void TestFitUsingClass(void)
{
    // Simple Test with 3 bins and two sources
    mystyle();
    TRandom3 * rd = new TRandom3(0);
    // Accept-Reject, using Root functionality. A-R w/ function
    
    // function: x^2 * exp(-x)
    
    int nIter = 1000000;
    double FactorGibbs = 0.1;
    double BurnInRatio = 0.2;
    int MaxLag=3000;
    
    TH1D * BinDefinition =  new TH1D("BinDefinition", "", 3, 0., 1.); // Empty histogram to clone
    
    Int_t nBins = 3;
    Double_t p1True = 0.5;
    Int_t n1 = 40;
    Double_t p2True = 0.7;
    Int_t n2 = 90;
    
    p1acc=0; // Global variables to check acceptance of acc-rej step in gibbs
    p1all=0;
    
    TH1D * Source1True = (TH1D*) BinDefinition->Clone("Source1True");
    Source1True->SetBinContent(1,0.25);
    Source1True->SetBinContent(2,1.5);
    Source1True->SetBinContent(3,0.25);
    TH1D * Source2True = (TH1D*) BinDefinition->Clone("Source2True");
    Source2True->SetBinContent(1,1./3.);
    Source2True->SetBinContent(2,1./3.);
    Source2True->SetBinContent(3,1./3.);
    
    Int_t nData = rd->Poisson(p1True*n1 + p2True*n2);
    TH1D * DataTrue = (TH1D*) BinDefinition->Clone("DataTrue");
    DataTrue->Add(Source1True, double(n1)*p1True/Source1True->Integral());
    DataTrue->Add(Source2True, double(n2)*p2True/Source2True->Integral());
    
    TH1D * DataHist = (TH1D*) BinDefinition->Clone("DataHist");
    DataHist->FillRandom(DataTrue, nData);
    TH1D * Template1 = (TH1D*) BinDefinition->Clone("Template1");
    Template1->FillRandom(Source1True, n1);
    TH1D * Template2 = (TH1D*) BinDefinition->Clone("Template2");
    Template2->FillRandom(Source2True, n2);
    
    A1iProposal = (TH1D*) BinDefinition->Clone("A1iProposal"); // Initialize temporary proposal points (should be in class later)
    A2iProposal = (TH1D*) BinDefinition->Clone("A2iProposal");
    pprop = new double[2];
    
    // State vectors for Metropolis
    TH1D * A1i = (TH1D*) Template1->Clone("A1i");
    TH1D * A2i = (TH1D*) Template2->Clone("A2i"); // Use templates for initial values
    double * p = new double[2];
    p[0] = 1.;
    p[1] = 1.;
    
    TH2D * pPosterior = new TH2D("pPosterior", "", 100., 0., 2., 100., 0., 2.); // pj marginal
    pPosterior->GetXaxis()->SetTitle("p_{1}");
    pPosterior->GetYaxis()->SetTitle("p_{2}");
    TH1D * p1Posterior = new TH1D("p1Posterior", "", 100., 0., 2.);
    p1Posterior->GetXaxis()->SetTitle("p_{1}");
    p1Posterior->GetYaxis()->SetTitle("normalized entries");
    TH1D * p2Posterior = new TH1D("p2Posterior", "", 100., 0., 2.);
    p2Posterior->GetXaxis()->SetTitle("p_{2}");
    p2Posterior->GetYaxis()->SetTitle("normalized entries");
    
    
    double * ValuesA11 = new double[nIter];
    double * ValuesA12 = new double[nIter];
    double * ValuesA13 = new double[nIter];
    double * ValuesA21 = new double[nIter];
    double * ValuesA22 = new double[nIter];
    double * ValuesA23 = new double[nIter];
    double * Valuesp1 = new double[nIter];
    double * Valuesp2 = new double[nIter]; // Metropolis Values for correlations
    
    double * Valuesp1Gibbs = new double[int(nIter*FactorGibbs)];
    double * Valuesp2Gibbs = new double[int(nIter*FactorGibbs)]; // Gibbs Values for correlations
    
    
    double TimeMetropolis, TimeGibbs;
    
    auto t1 = Clock::now();
    auto t2 = Clock::now();
    
    Int_t BurnInSteps = Int_t(nIter *BurnInRatio);
    cout << "Burn In with " << BurnInSteps << " steps\n0\%" << flush;
    for(int i=0;i<BurnInSteps;i++)
    {
        if((i*100)/BurnInSteps>((i-1)*100)/BurnInSteps) cout << "\r" << (i*100)/BurnInSteps << "\%  " << flush;
        DoMetropolisStep2Sources(Template1, Template2, DataHist, p, A1i, A2i, rd);
    }
    cout << "\r100\%" << endl;
    
    
    cout << "Metropolis MCMC with " << nIter << " steps\n0\%" << flush;
    Int_t accepted=0;
    t1 = Clock::now();
    for(int i=0;i<nIter;i++)
    {
        if((i*100)/nIter>((i-1)*100)/nIter) cout << "\r" << (i*100)/nIter << "\%  " << flush;
        if(DoMetropolisStep2Sources(Template1, Template2, DataHist, p, A1i, A2i, rd))
            accepted++;
        pPosterior->Fill(p[0],p[1]);
        p1Posterior->Fill(p[0]);
        p2Posterior->Fill(p[1]);
        ValuesA11[i] = A1i->GetBinContent(1);
        ValuesA12[i] = A1i->GetBinContent(2);
        ValuesA13[i] = A1i->GetBinContent(3);
        ValuesA21[i] = A2i->GetBinContent(1);
        ValuesA22[i] = A2i->GetBinContent(2);
        ValuesA23[i] = A2i->GetBinContent(3);
        Valuesp1[i] = p[0];
        Valuesp2[i] = p[1];
    
    }
    t2 = Clock::now();
    cout << "\r100\%" << endl;
    cout << "Acceptance ratio " << double(accepted)/double(nIter) << endl;
    TimeMetropolis = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "Time per iteration: " << TimeMetropolis/double(nIter) << " microseconds" << std::endl;
    
   
    // Now try with class
    TH1D ** MCHists = new TH1D*[2];
    MCHists[0]=Template1;
    MCHists[1]=Template2;
    AliMCMCTemplateFitter * fitClass = new AliMCMCTemplateFitter(2, MCHists, DataHist);
    fitClass->Fit(int(nIter*FactorGibbs));
    TH1D * p1GibbsClass = fitClass->ReturnpjMarginal(0);
    TH1D * p2GibbsClass = fitClass->ReturnpjMarginal(1);
    TH2D * pPosteriorGibbs = fitClass->ReturnpjMarginal(0,1);
    pPosteriorGibbs->SetAxisRange(0.,2.,"X");
    pPosteriorGibbs->SetAxisRange(0.,2.,"Y");
    p1GibbsClass->RebinX(10);
    p2GibbsClass->RebinX(10);
    
    
    cout << "Calculating correlations." << endl;
    
    TH1D * p1CorrelationGibbs = fitClass->GetpjAutocorrelation(0, 200);
    TH1D * p2CorrelationGibbs = fitClass->GetpjAutocorrelation(1, 200);
    
    TH1D * p1Correlation = new TH1D("p1Correlation", "", MaxLag, 0.5, 0.5+double(MaxLag));
    TH1D * p2Correlation = new TH1D("p2Correlation", "", MaxLag, 0.5, 0.5+double(MaxLag));
    
    double p1mean = p1Posterior->GetMean();
    double p2mean = p2Posterior->GetMean();
    double p1rms = p1Posterior->GetRMS();
    double p2rms = p2Posterior->GetRMS();
    double sump1sq, sump2sq;
    for(int lag=1;lag<=MaxLag;lag++)
    {
        sump1sq=0.;
        sump2sq=0.;
        for(int i=0;i<nIter-lag;i++)
        {
            sump1sq += (Valuesp1[i]-p1mean)*(Valuesp1[i+lag]-p1mean)/p1rms/p1rms;
            sump2sq += (Valuesp2[i]-p2mean)*(Valuesp2[i+lag]-p2mean)/p2rms/p2rms;
        }
        p1Correlation->SetBinContent(lag, sump1sq/double(nIter-lag));
        p2Correlation->SetBinContent(lag, sump2sq/double(nIter-lag));
    }
    
    
    
    /*TH1D * p1CorrelationGibbs = new TH1D("p1CorrelationGibbs", "", MaxLag, 0.5, 0.5+double(MaxLag));
    TH1D * p2CorrelationGibbs = new TH1D("p2CorrelationGibbs", "", MaxLag, 0.5, 0.5+double(MaxLag));
    
    double p1meanGibbs = p1Posterior->GetMean();
    double p2meanGibbs = p2Posterior->GetMean();
    double p1rmsGibbs = p1Posterior->GetRMS();
    double p2rmsGibbs = p2Posterior->GetRMS();
    for(int lag=1;lag<=MaxLag;lag++)
    {
        sump1sq=0.;
        sump2sq=0.;
        for(int i=0;i<int(nIter*FactorGibbs)-lag;i++)
        {
            sump1sq += (Valuesp1Gibbs[i]-p1meanGibbs)*(Valuesp1Gibbs[i+lag]-p1meanGibbs)/p1rmsGibbs/p1rmsGibbs;
            sump2sq += (Valuesp2Gibbs[i]-p2meanGibbs)*(Valuesp2Gibbs[i+lag]-p2meanGibbs)/p2rmsGibbs/p2rmsGibbs;
        }
        p1CorrelationGibbs->SetBinContent(lag, sump1sq/double(int(nIter*FactorGibbs)-lag));
        p2CorrelationGibbs->SetBinContent(lag, sump2sq/double(int(nIter*FactorGibbs)-lag));
    }*/
    cout << "Done." << endl;
    
    TH1D * BlankHisto = new TH1D("bg", "", 100, 0., 1.);
    BlankHisto->SetLineColor(kWhite);
    BlankHisto->SetAxisRange(0.1, DataHist->GetMaximum()*1.2,"Y");
    BlankHisto->GetXaxis()->SetTitle("x");
    BlankHisto->GetYaxis()->SetTitle("counts");
    //TimeAccRej = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    //std::cout << "Per iteration: " << TimeAccRej/double(nIter) << " microseconds" << std::endl;
    SetHistogram(DataHist, kBlack);
    SetHistogram(Template1, kGreen+1);
    SetHistogram(Template2, kRed+1);
    Template1->Scale(p1True);
    Template2->Scale(p2True);
    SetHistogram(p1Posterior, kBlack);
    SetHistogram(p2Posterior, kBlack);
    p1Posterior->Scale(1./p1Posterior->Integral("width"));
    p2Posterior->Scale(1./p2Posterior->Integral("width"));
    p1GibbsClass->Scale(1./p1GibbsClass->Integral("width"));
    p2GibbsClass->Scale(1./p2GibbsClass->Integral("width"));
    SetHistogram(p1GibbsClass, kGreen+1);
    SetHistogram(p2GibbsClass, kGreen+1);
    SetHistogram(p1Correlation, kBlack);
    SetHistogram(p2Correlation, kBlack);
    p2Correlation->SetLineStyle(2);
    p1Correlation->GetXaxis()->SetTitle("lag");
    p1Correlation->GetYaxis()->SetTitle("correlation coefficient");
    p1Correlation->SetAxisRange(0.01, 1.1, "Y");
    SetHistogram(p1CorrelationGibbs, kRed+1);
    SetHistogram(p2CorrelationGibbs, kRed+1);
    p2CorrelationGibbs->SetLineStyle(2);
    
    TMarker * Truth = new TMarker(p1True, p2True, 20);
    
    TLine * p1Truth = new TLine(p1True, 0., p1True, p1Posterior->GetMaximum()*0.1); p1Truth->SetLineWidth(2); p1Truth->SetLineColor(kBlue+1);
    TLine * p2Truth = new TLine(p2True, 0., p2True, p2Posterior->GetMaximum()*0.1); p2Truth->SetLineWidth(2); p2Truth->SetLineColor(kBlue+1);
    
    
    TCanvas * cv = new TCanvas("SampleTests", "", 1600, 1000);
    cv->Divide(3,2);
    cv->cd(1);
    BlankHisto->Draw();
    SetRect1DHist((TH1*)BlankHisto, (TPad*)gPad);
    DataHist->Draw("sameehist");
    Template1->Draw("sameehist");
    Template2->Draw("sameehist");
    TLegend * leg = plotLegend("right_top",  "", 0.9, 0.7, -0.05, 0., "",1);
    leg->AddEntry(DataHist, "Data", "l");
    leg->AddEntry(Template1, "Template1*p_{1,true}", "l");
    leg->AddEntry(Template2, "Template2*p_{2,true}", "l");
    leg->Draw("same");
    cv->cd(2);
    pPosterior->Draw("colz");
    TPaveText *pt = new TPaveText(0.3,0.7,0.7,0.8, "ndc");
    pt->AddText(Form("Metropolis, %d iterations", nIter));
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->Draw("same");
    SetSquare2DHist((TH2*)pPosterior, (TPad*)gPad);
    Truth->Draw("same");
    cv->cd(3);
    pPosteriorGibbs->Draw("colz");
    TPaveText *ptGibbs = new TPaveText(0.3,0.7,0.7,0.8, "ndc");
    ptGibbs->AddText(Form("Gibbs, %d iterations", int(nIter*FactorGibbs)));
    ptGibbs->SetFillStyle(0);
    ptGibbs->SetBorderSize(0);
    ptGibbs->Draw("same");
    SetSquare2DHist((TH2*)pPosteriorGibbs, (TPad*)gPad);
    Truth->Draw("same");
    //BlankHisto->Draw("p");
    
    cv->cd(4);
    p1Posterior->Draw("hist");
    SetRect1DHist((TH1*)p1Posterior, (TPad*)gPad);
    p1GibbsClass->Draw("samehist");
    p1Truth->Draw("same");
    TLegend * legp1 = plotLegend("right_top",  "", 0.9, 0.5, 0., 0., "p_{1} marginal",1);
    legp1->AddEntry(p1Posterior, "Metropolis", "l");
    legp1->AddEntry(p1GibbsClass, "Gibbs", "l");
    legp1->AddEntry(p1Truth, "True Value", "l");
    legp1->Draw("same");
    cv->cd(5);
    p2Posterior->Draw("hist");
    SetRect1DHist((TH1*)p2Posterior, (TPad*)gPad);
    p2GibbsClass->Draw("samehist");
    p2Truth->Draw("same");
    TLegend * legp2 = plotLegend("right_top",  "", 0.9, 0.5, 0., 0., "p_{2} marginal",1);
    legp2->AddEntry(p2Posterior, "Metropolis", "l");
    legp2->AddEntry(p2GibbsClass, "Gibbs", "l");
    legp2->AddEntry(p2Truth, "True Value", "l");
    legp2->Draw("same");
    cv->cd(6);
    p1Correlation->Draw();
    SetRect1DHist((TH1*)p1Correlation, (TPad*)gPad);
    p2Correlation->Draw("same");
    p1CorrelationGibbs->Draw("same");
    p2CorrelationGibbs->Draw("same");
    gPad->SetLogx();
    gPad->SetLogy();
    TLegend * legcorr = plotLegend("right_top",  "", 0.9, 0.5, 0., 0., "Autocorrelation",1);
    legcorr->AddEntry(p1Correlation, "Metropolis p_{1}", "l");
    legcorr->AddEntry(p2Correlation, "Metropolis p_{2}", "l");
    legcorr->AddEntry(p1CorrelationGibbs, "Gibbs p_{1}", "l");
    legcorr->AddEntry(p2CorrelationGibbs, "Gibbs p_{2}", "l");
    legcorr->Draw("same");
}

