#include "MartinsStyle.h"
#include "AliMCMCTemplateFitter.h"
#include "AliMCLogLFitter.h"
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

void SetHistogram(TH1D * in, int color)
{
    in->SetLineWidth(2);
    in->SetLineColor(color);
    in->SetMarkerColor(color);
}

void MakeTruthMonotonous(TH1D* in)
{
    int n=in->GetNbinsX();
    for(int i=n/2-2;i>0;i--)
        if(in->GetBinContent(i) > in->GetBinContent(i+1)) in->SetBinContent(i, in->GetBinContent(i+1)*0.95);
    for(int i=n/2+2;i<=n;i++)
        if(in->GetBinContent(i) > in->GetBinContent(i-1)) in->SetBinContent(i, in->GetBinContent(i-1)*0.95); 
}

void SetHistoFct(TH1D * h, TF1 * f)
{
    int n = h->GetNbinsX();
    for(int i=1;i<=n;i++)
        h->SetBinContent(i, f->Eval(h->GetXaxis()->GetBinCenter(i)));
}

void TestFitUsingClass(void)
{
    // Simple Test with 3 bins and two sources
    mystyle();
    TRandom3 * rd = new TRandom3(7);
    gRandom->SetSeed(3);
    cout << "RNG seed: " << gRandom->GetSeed() << endl;
    // Accept-Reject, using Root functionality. A-R w/ function
    
    // function: x^2 * exp(-x)
    
    int nIter = 3000;
    Bool_t MakeAjiAutocorr = kTRUE;
    int nBins = 100;
    
    Double_t p0True = 1.3;
    Int_t n0 = 200;
    Double_t p1True = 1.4;
    Int_t n1 = 200;
    Double_t p2True = 1.5;
    Int_t n2 = 200;
    
    TF1 * fct0 = new TF1("fct0", "TMath::Exp(-TMath::Abs(x/0.2))", -1., 1.);
    TF1 * fct1 = new TF1("fct1", "TMath::Exp(-TMath::Abs(x/0.5))", -1., 1.);
    TF1 * fct2 = new TF1("fct2", "TMath::Exp(-TMath::Abs((x-0.3)/0.2))", -1., 1.);
    
    TH1D * BinDef = new TH1D("BinDef", "", nBins, -1., 1.);
    
    TH1D * H0 =  (TH1D*) BinDef->Clone("H0");
    TH1D * H1 =  (TH1D*) BinDef->Clone("H1");
    TH1D * H2 =  (TH1D*) BinDef->Clone("H2");
    
    TH1D * Template0 = (TH1D*) BinDef->Clone("Template0");
    TH1D * Template1 = (TH1D*) BinDef->Clone("Template1");
    TH1D * Template2 = (TH1D*) BinDef->Clone("Template2");
    
    SetHistoFct(H0, fct0);
    SetHistoFct(H1, fct1);
    SetHistoFct(H2, fct2);
    
    Int_t nData = rd->Poisson(p0True*double(n0) + p1True*double(n1) + p2True*double(n2));
    TH1D * DataTrue = (TH1D*) BinDef->Clone("DataTrue");
    DataTrue->Add(H0, double(n0)*p0True/H0->Integral());
    DataTrue->Add(H1, double(n0)*p1True/H1->Integral());
    DataTrue->Add(H2, double(n0)*p2True/H2->Integral());
    
    TH1D * DataHist = (TH1D*) BinDef->Clone("DataHist");
    DataHist->FillRandom(DataTrue, nData);
    Template0->FillRandom(H0, n0);
    Template1->FillRandom(H1, n1);
    Template2->FillRandom(H2, n2);
    cout << "Maximum template1: " << Template0->GetMaximum() << " template2: " << Template1->GetMaximum() << " template3: " << Template2->GetMaximum() <<  " data: " << DataHist->GetMaximum() <<  endl;
    cout << "Largest bin : " << H0->GetMaximumBin() << endl;
    
    
    // Now try with class
    TH1D ** MCHists = new TH1D*[3];
    MCHists[0]=Template0;
    MCHists[1]=Template1;
    MCHists[2]=Template2;
    AliMCMCTemplateFitter * fitClass = new AliMCMCTemplateFitter(3, MCHists, DataHist);
    if(MakeAjiAutocorr) fitClass->SetPrepareAjiAutocorrelation();
    fitClass->SetBurnInRatio(0.5);
    //fitClass->SetMonotonyPrior(0, 49, 50);
    //fitClass->SetMonotonyPrior(1, 49, 50);
    //fitClass->SetMonotonyPrior(2, 64, 65);
    //fitClass->SetFitRange(-0.1, 0.1); // include limits also in monotony prior
    fitClass->Fit(nIter);
    TH1D * p1GibbsClass = fitClass->ReturnpjMarginal(0);
    TH1D * p2GibbsClass = fitClass->ReturnpjMarginal(1);
    TH1D * p3GibbsClass = fitClass->ReturnpjMarginal(2);
    TH2D * pPosteriorGibbs = fitClass->ReturnpjMarginal(0,1);
    TH2D * pPosteriorGibbs12 = fitClass->ReturnpjMarginal(1,2);
    p1GibbsClass->RebinX(5);
    p2GibbsClass->RebinX(5);
    p3GibbsClass->RebinX(5);
    
    
    cout << "Calculating correlations." << endl;
    
    int MaxCorr = 500;
    
    TH1D * p1CorrelationGibbs = fitClass->GetpjAutocorrelation(0, MaxCorr);
    TH1D * p2CorrelationGibbs = fitClass->GetpjAutocorrelation(1, MaxCorr);
    TH1D * p3CorrelationGibbs = fitClass->GetpjAutocorrelation(2, MaxCorr);
    
    TH1D ** AjiCorrelations = new TH1D*[3*nBins];
    
    TH2D * A0iMarginal = fitClass->ReturnAjiMarginals(0);
    TH2D * A1iMarginal = fitClass->ReturnAjiMarginals(1);
    TH2D * A2iMarginal = fitClass->ReturnAjiMarginals(2);
    if(MakeAjiAutocorr)
    for(int j = 0;j<3;j++)
        for(int i=0;i<nBins;i++)
    
    SetHistogram(H0, kGreen+1);
    SetHistogram(H1, kRed+1);
    SetHistogram(H2, kGray+1);
    H0->Scale(1./H0->Integral());
    H1->Scale(1./H1->Integral());
    H2->Scale(1./H2->Integral());
    double low= H0->GetXaxis()->GetBinLowEdge(1);
    double high = H0->GetXaxis()->GetBinLowEdge(nBins+1);
    cout << "3 * " << nBins << " bins. Leading to " << 3*nBins +3 << " unknown free parameters" << endl;
    
    for(int j = 0;j<3;j++)
        for(int i=0;i<nBins;i++)
        {
            AjiCorrelations[i*3+j] = fitClass->GetAjiAutocorrelation(j, i, MaxCorr);
            AjiCorrelations[i*3+j]->SetLineColor(kBlack);
        }
    
    
    cout << "Done." << endl;
    
    cout << "Now checking with the ML fit." << endl;
    AliMCLogLFitter * MLFitter = new AliMCLogLFitter(3, (TH1**)MCHists, (TH1*)DataHist);
    //Fitter->SetUseChiSq();
    MLFitter->IfNoMCParticleIsProbablyA(1);
    //Fitter->SetFitRange(-0.10, 0.10);
    MLFitter->Fit();
    double MLpar0 = MLFitter->ReturnParameter(0);
    double MLpar1 = MLFitter->ReturnParameter(1);
    double MLpar2 = MLFitter->ReturnParameter(2);
    TMarker * MLMarker01 = new TMarker(MLpar0, MLpar1, 20); MLMarker01->SetMarkerColor(kRed+2);
    TMarker * MLMarker12 = new TMarker(MLpar1, MLpar2, 20); MLMarker12->SetMarkerColor(kRed+2);
    cout << "end of ML fit" << endl;
    
    SetHistogram(DataHist, kBlack);
    SetHistogram(Template0, kGreen+1);
    SetHistogram(Template1, kRed+1);
    TH1D * H0InputCopy = (TH1D*) Template0->Clone("H0InputCopy");
    TH1D * H1InputCopy = (TH1D*) Template1->Clone("H1InputCopy");
    TH1D * ConversionInputCopy = (TH1D*) Template2->Clone("ConversionInputCopy");
    Template0->Scale(p1GibbsClass->GetMean());
    Template1->Scale(p2GibbsClass->GetMean());
    Template2->Scale(p3GibbsClass->GetMean());
    TH1D * Sum = (TH1D*) BinDef->Clone("Sum");
    Sum->Add(Template0);
    Sum->Add(Template1);
    Sum->Add(Template2);
    Sum->SetLineWidth(2);Sum->SetLineColor(kOrange+1);
    p1GibbsClass->Scale(1./p1GibbsClass->Integral("width"));
    p2GibbsClass->Scale(1./p2GibbsClass->Integral("width"));
    p3GibbsClass->Scale(1./p3GibbsClass->Integral("width"));
    SetHistogram(p1GibbsClass, kGreen+1);
    SetHistogram(p2GibbsClass, kRed+1);
    SetHistogram(p3GibbsClass, kGray+1);
    SetHistogram(p1CorrelationGibbs, kGreen+1);
    SetHistogram(p2CorrelationGibbs, kRed+1);
    SetHistogram(p3CorrelationGibbs, kGray+1);
    p2CorrelationGibbs->SetLineStyle(2);
    p1CorrelationGibbs->SetAxisRange(0.07, 1.1, "Y");
    
    TMarker * Truth = new TMarker(p0True, p1True, 20);
    TMarker * Truth12 = new TMarker(p1True, p2True, 20);
    
    TH2D * cont12 = fitClass->GetCredibleArea(1, 2, 0.95);
    TH2D * cont01 = fitClass->GetCredibleArea(0, 1, 0.95);
    
    double MarginalMaximum = TMath::Max(p1GibbsClass->GetMaximum(), p2GibbsClass->GetMaximum());
    p2GibbsClass->SetAxisRange(0., MarginalMaximum*1.2, "Y");
    
    TLine * p1Truth = new TLine(p0True, 0., p0True, MarginalMaximum*0.1); p1Truth->SetLineWidth(2); p1Truth->SetLineColor(kGreen+1);
    TLine * p2Truth = new TLine(p1True, 0., p1True, MarginalMaximum*0.1); p2Truth->SetLineWidth(2); p2Truth->SetLineColor(kRed+1);
    TLine * p3Truth = new TLine(p2True, 0., p2True, MarginalMaximum*0.1); p3Truth->SetLineWidth(2); p3Truth->SetLineColor(kGray+1);
    
    TH1D * BlankHisto = new TH1D("bg", "", 100, -1., 1.);
    BlankHisto->SetLineColor(kWhite);
    BlankHisto->SetAxisRange(0.1, DataHist->GetMaximum()*1.2,"Y");
    BlankHisto->GetXaxis()->SetTitle("x");
    BlankHisto->GetYaxis()->SetTitle("counts");
        
    
    
    TCanvas * cv = new TCanvas("IPFitExample", "", 1600, 1000);
    cv->Divide(3,2);
    cv->cd(1);
    BlankHisto->Draw();
    SetRect1DHist((TH1*)BlankHisto, (TPad*)gPad);
    //A1iMarginal->Draw("samecolz");
    Sum->Draw("sameehist");
    DataHist->Draw("sameehist");
    Template0->Draw("sameehist");
    Template1->Draw("sameehist");
    Template2->Draw("sameehist");
    gPad->SetLogy();
    TLegend * leg = plotLegend("right_top",  "", 0.9, 0.7, -0.05, 0., "",1);
    leg->AddEntry(DataHist, "Data", "l");
    leg->AddEntry(Template0, "Template0*p_{0}", "l");
    leg->AddEntry(Template1, "Template1*p_{1}", "l");
    leg->AddEntry(Template2, "Template2*p_{2}", "l");
    leg->AddEntry(Sum, "Sum of scaled dist", "l");
    leg->Draw("same");
    cv->cd(2);
    pPosteriorGibbs->Draw("colz");
    SetSquare2DHist((TH2*)pPosteriorGibbs, (TPad*)gPad);
    cont01->Draw("samecont2");
    Truth->Draw("same");
    MLMarker01->Draw("same");
    TLegend * leg01 = plotLegend("right_top",  "", 1.2, 0.7, -0.2, 0., Form("Gibbs, %d iterations", nIter),1);
    leg01->AddEntry(Truth, "True value", "p");
    leg01->AddEntry(cont01, "95% credible area", "l");
    leg01->AddEntry(MLMarker12, "ML", "p");
    leg01->Draw("same");
    cv->cd(3);
    pPosteriorGibbs12->Draw("colz");
    SetSquare2DHist((TH2*)pPosteriorGibbs12, (TPad*)gPad);
    Truth12->Draw("same");
    MLMarker12->Draw("same");
    cont12->Draw("samecont2");
    cv->cd(4);
    H0->Draw("hist");
    SetRect1DHist((TH1*)H0, (TPad*)gPad);
    H1->Draw("SAMEhist");
    H2->Draw("samehist");
    gPad->SetLogy();
    cv->cd(5);
    p2GibbsClass->SetAxisRange(0., TMath::Max(TMath::Max(p1GibbsClass->GetMaximum(), p2GibbsClass->GetMaximum()), p3GibbsClass->GetMaximum())*1.2, "Y");
    p2GibbsClass->Draw("hist");
    p1GibbsClass->Draw("samehist");
    p3GibbsClass->Draw("samehist");
    SetRect1DHist((TH1*)p2GibbsClass, (TPad*)gPad);
    p2Truth->Draw("same");
    p1Truth->Draw("same");
    p3Truth->Draw("same");
    TLegend * legp2 = plotLegend("right_top",  "", 1.1, 0.8, -0.1, 0., "marginals with true values",1);
    legp2->AddEntry(p1GibbsClass, "p_{0}", "l");
    legp2->AddEntry(p2GibbsClass, "p_{1}", "l");
    legp2->AddEntry(p3GibbsClass, "p_{2}", "l");
    legp2->Draw("same");
    cv->cd(6);
    p1CorrelationGibbs->Draw("");
    SetRect1DHist((TH1*)p1CorrelationGibbs, (TPad*)gPad);
    p2CorrelationGibbs->Draw("same");
    p3CorrelationGibbs->Draw("same");
    if(MakeAjiAutocorr)
    for(int j = 0;j<3;j++)
        for(int i=0;i<nBins;i++)
            AjiCorrelations[i*3+j]->Draw("same");;
    gPad->SetLogx();
    gPad->SetLogy();
    TLegend * legcorr = plotLegend("right_top",  "", 0.9, 0.6, 0., 0., "Autocorrelation",1);
    legcorr->AddEntry(p1CorrelationGibbs, "p_{0}", "l");
    legcorr->AddEntry(p2CorrelationGibbs, "p_{1}", "l");
    legcorr->AddEntry(p3CorrelationGibbs, "p_{2}", "l");
    //legcorr->AddEntry(AjiCorrelations[0], "Gibbs A_{ji}s", "l");
    legcorr->Draw("same");
    
    
    TCanvas * cv2 = new TCanvas("H1AjiMarginals", "", 1400, 1000);
    cv2->Divide(2,2);
    cv2->cd(1);
    A0iMarginal->Draw("colz");
    H0InputCopy->Draw("histesame");
    gPad->SetLogy();
    gPad->SetLogz();
    cv2->cd(2);
    A1iMarginal->Draw("colz");
    H1InputCopy->Draw("histesame");
    gPad->SetLogy();
    gPad->SetLogz();
    cv2->cd(3);
    A2iMarginal->Draw("colz");
    ConversionInputCopy->Draw("histesame");
    gPad->SetLogy();
    gPad->SetLogz();
    cv2->cd(4);
    
}

