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


void TestFitUsingClass(void)
{
    // Simple Test with 3 bins and two sources
    mystyle();
    TRandom3 * rd = new TRandom3(0);
    // Accept-Reject, using Root functionality. A-R w/ function
    
    // function: x^2 * exp(-x)
    
    int nIter = 10000;
    
    double low= -1.;
    double high = 1.;
    
    Int_t nBins = 50;
    Double_t p1True = 0.45;
    Int_t n1 = 1290;
    Double_t p2True = 0.26;
    Int_t n2 = 2100;
    
    TH1D * BinDefinition =  new TH1D("BinDefinition", "", nBins, low, high); // Empty histogram to clone
    
    TF1 * f1 = new TF1("f1", "TMath::Exp(-TMath::Abs(3.*x))", low, high);
    //TF1 * f2 = new TF1("f1", "TMath::Exp(-x^2/(2.*1.5*1.5))", low, high);
    TF1 * f2 = new TF1("f2", "TMath::Exp(-TMath::Abs(0.9*x))", low, high);
    
    
    TH1D * Source1True = (TH1D*) BinDefinition->Clone("Source1True");
    for(int i=1;i<=nBins;i++)
        Source1True->SetBinContent(i, f1->Eval(Source1True->GetXaxis()->GetBinCenter(i)));
    TH1D * Source2True = (TH1D*) BinDefinition->Clone("Source2True");
    for(int i=1;i<=nBins;i++)
        Source2True->SetBinContent(i, f2->Eval(Source2True->GetXaxis()->GetBinCenter(i)));
    
    Int_t nData = rd->Poisson(p1True*double(n1) + p2True*double(n2));
    TH1D * DataTrue = (TH1D*) BinDefinition->Clone("DataTrue");
    DataTrue->Add(Source1True, double(n1)*p1True/Source1True->Integral());
    DataTrue->Add(Source2True, double(n2)*p2True/Source2True->Integral());
    
    TH1D * DataHist = (TH1D*) BinDefinition->Clone("DataHist");
    DataHist->FillRandom(DataTrue, nData);
    TH1D * Template1 = (TH1D*) BinDefinition->Clone("Template1");
    Template1->FillRandom(Source1True, n1);
    TH1D * Template2 = (TH1D*) BinDefinition->Clone("Template2");
    Template2->FillRandom(Source2True, n2);
    cout << "Maximum template1: " << Template1->GetMaximum() << " template2: " << Template2->GetMaximum() << " data: " << DataHist->GetMaximum() << endl;
    
    
    // Now try with class
    TH1D ** MCHists = new TH1D*[2];
    MCHists[0]=Template1;
    MCHists[1]=Template2;
    AliMCMCTemplateFitter * fitClass = new AliMCMCTemplateFitter(2, MCHists, DataHist);
    fitClass->Fit(nIter);
    TH1D * p1GibbsClass = fitClass->ReturnpjMarginal(0);
    TH1D * p2GibbsClass = fitClass->ReturnpjMarginal(1);
    TH2D * pPosteriorGibbs = fitClass->ReturnpjMarginal(0,1);
    //pPosteriorGibbs->SetAxisRange(0.,2.,"X");
    //pPosteriorGibbs->SetAxisRange(0.,2.,"Y");
    p1GibbsClass->RebinX(10);
    p2GibbsClass->RebinX(10);
    
    
    cout << "Calculating correlations." << endl;
    
    TH1D * p1CorrelationGibbs = fitClass->GetpjAutocorrelation(0, 200);
    TH1D * p2CorrelationGibbs = fitClass->GetpjAutocorrelation(1, 200);
    
    cout << "Done." << endl;
    
    TH1D * BlankHisto = new TH1D("bg", "", 100, low, high);
    BlankHisto->SetLineColor(kWhite);
    BlankHisto->SetAxisRange(0., DataHist->GetMaximum()*1.2,"Y");
    BlankHisto->GetXaxis()->SetTitle("x");
    BlankHisto->GetYaxis()->SetTitle("counts");
    SetHistogram(DataHist, kBlack);
    SetHistogram(Template1, kGreen+1);
    SetHistogram(Template2, kRed+1);
    Template1->Scale(p1True);
    Template2->Scale(p2True);
    TH1D * Sum = (TH1D*) BinDefinition->Clone("Sum");
    Sum->Add(Template1);
    Sum->Add(Template2);
    Sum->SetLineWidth(2);Sum->SetLineColor(kBlue+1);
    p1GibbsClass->Scale(1./p1GibbsClass->Integral("width"));
    p2GibbsClass->Scale(1./p2GibbsClass->Integral("width"));
    SetHistogram(p1GibbsClass, kGreen+1);
    SetHistogram(p2GibbsClass, kGreen+1);
    SetHistogram(p1CorrelationGibbs, kRed+1);
    SetHistogram(p2CorrelationGibbs, kRed+1);
    p2CorrelationGibbs->SetLineStyle(2);
    p1CorrelationGibbs->SetAxisRange(0.01, 1.1, "Y");
    
    TMarker * Truth = new TMarker(p1True, p2True, 20);
    
    TLine * p1Truth = new TLine(p1True, 0., p1True, p1GibbsClass->GetMaximum()*0.1); p1Truth->SetLineWidth(2); p1Truth->SetLineColor(kBlue+1);
    TLine * p2Truth = new TLine(p2True, 0., p2True, p2GibbsClass->GetMaximum()*0.1); p2Truth->SetLineWidth(2); p2Truth->SetLineColor(kBlue+1);
    
    
    TCanvas * cv = new TCanvas("SampleTests", "", 1600, 1000);
    cv->Divide(3,2);
    cv->cd(1);
    BlankHisto->Draw();
    SetRect1DHist((TH1*)BlankHisto, (TPad*)gPad);
    Sum->Draw("sameehist");
    DataHist->Draw("sameehist");
    Template1->Draw("sameehist");
    Template2->Draw("sameehist");
    TLegend * leg = plotLegend("right_top",  "", 0.9, 0.7, -0.05, 0., "",1);
    leg->AddEntry(DataHist, "Data", "l");
    leg->AddEntry(Template1, "Template1*p_{0,true}", "l");
    leg->AddEntry(Template2, "Template2*p_{1,true}", "l");
    leg->AddEntry(Sum, "Sum of both", "l");
    leg->Draw("same");
    cv->cd(2);
    pPosteriorGibbs->Draw("colz");
    TPaveText *ptGibbs = new TPaveText(0.3,0.7,0.7,0.8, "ndc");
    ptGibbs->AddText(Form("Gibbs, %d iterations", nIter));
    ptGibbs->SetFillStyle(0);
    ptGibbs->SetBorderSize(0);
    ptGibbs->Draw("same");
    SetSquare2DHist((TH2*)pPosteriorGibbs, (TPad*)gPad);
    Truth->Draw("same");
    cv->cd(3);
    f1->SetLineColor(kGreen+1);
    f1->Draw();
    Source1True->Draw("samehist");
    Source2True->Draw("samehist");
    f2->Draw("same");
    cv->cd(4);
    p1GibbsClass->Draw("hist");
    SetRect1DHist((TH1*)p1GibbsClass, (TPad*)gPad);
    p1Truth->Draw("same");
    TLegend * legp1 = plotLegend("right_top",  "", 0.9, 0.5, 0., 0., "p_{0} marginal",1);
    legp1->AddEntry(p1GibbsClass, "Gibbs", "l");
    legp1->AddEntry(p1Truth, "True Value", "l");
    legp1->Draw("same");
    cv->cd(5);
    p2GibbsClass->Draw("hist");
    SetRect1DHist((TH1*)p2GibbsClass, (TPad*)gPad);
    p2Truth->Draw("same");
    TLegend * legp2 = plotLegend("right_top",  "", 0.9, 0.5, 0., 0., "p_{1} marginal",1);
    legp2->AddEntry(p2GibbsClass, "Gibbs", "l");
    legp2->AddEntry(p2Truth, "True Value", "l");
    legp2->Draw("same");
    cv->cd(6);
    p1CorrelationGibbs->Draw("");
    SetRect1DHist((TH1*)p1CorrelationGibbs, (TPad*)gPad);
    p2CorrelationGibbs->Draw("same");
    gPad->SetLogx();
    gPad->SetLogy();
    TLegend * legcorr = plotLegend("right_top",  "", 0.9, 0.5, 0., 0., "Autocorrelation",1);
    legcorr->AddEntry(p1CorrelationGibbs, "Gibbs p_{0}", "l");
    legcorr->AddEntry(p2CorrelationGibbs, "Gibbs p_{1}", "l");
    legcorr->Draw("same");
}

