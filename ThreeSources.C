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
    
    int nIter = 50000;
    
    double low= -1.;
    double high = 1.;
    
    Int_t nBins = 15;
    Double_t p1True = 0.62;
    Int_t n1 = 250;
    Double_t p2True = 0.85;
    Int_t n2 = 310;
    Double_t p3True = 0.51;
    Int_t n3 = 294;
    
    TH1D * BinDefinition =  new TH1D("BinDefinition", "", nBins, low, high); // Empty histogram to clone
    
    TF1 * f1 = new TF1("f1", "TMath::Exp(-TMath::Abs(3.*x))", low, high);
    //TF1 * f2 = new TF1("f1", "TMath::Exp(-x^2/(2.*1.5*1.5))", low, high);
    TF1 * f2 = new TF1("f2", "TMath::Exp(-TMath::Abs(0.8*x))", low, high);
    TF1 * f3 = new TF1("f3", "0.5- 0.1* x", low, high);
    
    
    TH1D * Source1True = (TH1D*) BinDefinition->Clone("Source1True");
    for(int i=1;i<=nBins;i++)
        Source1True->SetBinContent(i, f1->Eval(Source1True->GetXaxis()->GetBinCenter(i)));
    TH1D * Source2True = (TH1D*) BinDefinition->Clone("Source2True");
    for(int i=1;i<=nBins;i++)
        Source2True->SetBinContent(i, f2->Eval(Source2True->GetXaxis()->GetBinCenter(i)));
    TH1D * Source3True = (TH1D*) BinDefinition->Clone("Source3True");
    for(int i=1;i<=nBins;i++)
        Source3True->SetBinContent(i, f3->Eval(Source3True->GetXaxis()->GetBinCenter(i)));
    
    Int_t nData = rd->Poisson(p1True*double(n1) + p2True*double(n2) + p3True*double(n3));
    TH1D * DataTrue = (TH1D*) BinDefinition->Clone("DataTrue");
    DataTrue->Add(Source1True, double(n1)*p1True/Source1True->Integral());
    DataTrue->Add(Source2True, double(n2)*p2True/Source2True->Integral());
    DataTrue->Add(Source3True, double(n3)*p3True/Source3True->Integral());
    
    TH1D * DataHist = (TH1D*) BinDefinition->Clone("DataHist");
    DataHist->FillRandom(DataTrue, nData);
    TH1D * Template1 = (TH1D*) BinDefinition->Clone("Template1");
    Template1->FillRandom(Source1True, n1);
    TH1D * Template2 = (TH1D*) BinDefinition->Clone("Template2");
    Template2->FillRandom(Source2True, n2);
    TH1D * Template3 = (TH1D*) BinDefinition->Clone("Template3");
    Template3->FillRandom(Source3True, n3);
    cout << "Maximum template1: " << Template1->GetMaximum() << " template2: " << Template2->GetMaximum() << " template3: " << Template3->GetMaximum() << " data: " << DataHist->GetMaximum() << endl;
    
    
    // Now try with class
    TH1D ** MCHists = new TH1D*[3];
    MCHists[0]=Template1;
    MCHists[1]=Template2;
    MCHists[2]=Template3;
    AliMCMCTemplateFitter * fitClass = new AliMCMCTemplateFitter(3, MCHists, DataHist);
    fitClass->Fit(nIter);
    TH1D * p1GibbsClass = fitClass->ReturnpjMarginal(0);
    TH1D * p2GibbsClass = fitClass->ReturnpjMarginal(1);
    TH1D * p3GibbsClass = fitClass->ReturnpjMarginal(2);
    TH2D * pPosteriorGibbs = fitClass->ReturnpjMarginal(0,1);
    TH2D * pPosteriorGibbs12 = fitClass->ReturnpjMarginal(1,2);
    //pPosteriorGibbs->SetAxisRange(0.,2.,"X");
    //pPosteriorGibbs->SetAxisRange(0.,2.,"Y");
    p1GibbsClass->RebinX(10);
    p2GibbsClass->RebinX(10);
    p3GibbsClass->RebinX(10);
    
    
    cout << "Calculating correlations." << endl;
    
    TH1D * p1CorrelationGibbs = fitClass->GetpjAutocorrelation(0, 300);
    TH1D * p2CorrelationGibbs = fitClass->GetpjAutocorrelation(1, 300);
    TH1D * p3CorrelationGibbs = fitClass->GetpjAutocorrelation(2, 300);
    
    TH1D ** AjiCorrelations = new TH1D*[3*nBins];
    
    TH2D * A0iMarginal = fitClass->ReturnAjiMarginals(0);
    //for(int j = 0;j<3;j++)
    //    for(int i=0;i<nBins;i++)
    //    {
    //        AjiCorrelations[i*3+j] = fitClass->GetAjiAutocorrelation(j, i, 300);
    //        AjiCorrelations[i*3+j]->SetLineColor(kBlack);
    //    }
    
    
    cout << "Done." << endl;
    
    TH1D * BlankHisto = new TH1D("bg", "", 100, low, high);
    BlankHisto->SetLineColor(kWhite);
    BlankHisto->SetAxisRange(0., DataHist->GetMaximum()*1.2,"Y");
    BlankHisto->GetXaxis()->SetTitle("x");
    BlankHisto->GetYaxis()->SetTitle("counts");
    SetHistogram(DataHist, kBlack);
    SetHistogram(Template1, kGreen+1);
    SetHistogram(Template2, kRed+1);
    SetHistogram(Template3, kGray+1);
    //Template1->Scale(p1True);
    Template2->Scale(p2True);
    Template3->Scale(p3True);
    TH1D * Sum = (TH1D*) BinDefinition->Clone("Sum");
    Sum->Add(Template1, p1True);
    Sum->Add(Template2);
    Sum->Add(Template3);
    Sum->SetLineWidth(2);Sum->SetLineColor(kBlue+1);
    p1GibbsClass->Scale(1./p1GibbsClass->Integral("width"));
    p2GibbsClass->Scale(1./p2GibbsClass->Integral("width"));
    p3GibbsClass->Scale(1./p3GibbsClass->Integral("width"));
    SetHistogram(p1GibbsClass, kGreen+1);
    SetHistogram(p2GibbsClass, kRed+1);
    SetHistogram(p3GibbsClass, kGray+1);
    SetHistogram(p1CorrelationGibbs, kRed+1);
    SetHistogram(p2CorrelationGibbs, kGreen+1);
    SetHistogram(p3CorrelationGibbs, kGray+1);
    p2CorrelationGibbs->SetLineStyle(2);
    p1CorrelationGibbs->SetAxisRange(0.01, 1.1, "Y");
    
    TMarker * Truth = new TMarker(p1True, p2True, 20);
    TMarker * Truth12 = new TMarker(p2True, p3True, 20);
    
    TH2D * cont12 = fitClass->GetCredibleArea(1, 2, 0.9);
    TH2D * cont01 = fitClass->GetCredibleArea(0, 1, 0.9);
    
    double MarginalMaximum = TMath::Max(TMath::Max(p1GibbsClass->GetMaximum(), p2GibbsClass->GetMaximum()), p3GibbsClass->GetMaximum());
    p2GibbsClass->SetAxisRange(0., MarginalMaximum*1.2, "Y");
    
    TLine * p1Truth = new TLine(p1True, 0., p1True, MarginalMaximum*0.1); p1Truth->SetLineWidth(2); p1Truth->SetLineColor(kGreen+1);
    TLine * p2Truth = new TLine(p2True, 0., p2True, MarginalMaximum*0.1); p2Truth->SetLineWidth(2); p2Truth->SetLineColor(kRed+1);
    TLine * p3Truth = new TLine(p3True, 0., p3True, MarginalMaximum*0.1); p3Truth->SetLineWidth(2); p3Truth->SetLineColor(kGray+1);
    
    
    TCanvas * cv = new TCanvas("ThreeSourceExample", "", 1600, 1000);
    cv->Divide(3,2);
    cv->cd(1);
    BlankHisto->Draw();
    SetRect1DHist((TH1*)BlankHisto, (TPad*)gPad);
    A0iMarginal->Draw("samecolz");
    Sum->Draw("sameehist");
    DataHist->Draw("sameehist");
    Template1->Draw("sameehist");
    Template2->Draw("sameehist");
    Template3->Draw("sameehist");
    TLegend * leg = plotLegend("right_top",  "", 0.9, 0.7, -0.05, 0., "",1);
    leg->AddEntry(DataHist, "Data", "l");
    leg->AddEntry(Template1, "Template0*p_{0,true}", "l");
    leg->AddEntry(Template2, "Template1*p_{1,true}", "l");
    leg->AddEntry(Template3, "Template2*p_{2,true}", "l");
    leg->AddEntry(Sum, "Sum of both", "l");
    leg->Draw("same");
    cv->cd(2);
    pPosteriorGibbs->Draw("colz");
    SetSquare2DHist((TH2*)pPosteriorGibbs, (TPad*)gPad);
    cont01->Draw("samecont2");
    Truth->Draw("same");
    TLegend * leg01 = plotLegend("right_top",  "", 1.2, 0.7, -0.2, 0., Form("Gibbs, %d iterations", nIter),1);
    leg01->AddEntry(Truth, "True value", "p");
    leg01->AddEntry(cont01, "95% credible area", "l");
    leg01->Draw("same");
    cv->cd(3);
    pPosteriorGibbs12->Draw("colz");
    SetSquare2DHist((TH2*)pPosteriorGibbs12, (TPad*)gPad);
    Truth12->Draw("same");
    cont12->Draw("samecont2");
    cv->cd(4);
    Source1True->Draw("hist");
    SetRect1DHist((TH1*)Source1True, (TPad*)gPad);
    f1->SetLineColor(kGreen+1);
    f1->Draw("SAME");
    Source2True->Draw("samehist");
    Source3True->Draw("samehist");
    f2->Draw("same");
    f3->SetLineColor(kGray+1);
    f3->Draw("same");
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
    p3CorrelationGibbs->Draw("same");for(int j = 0;j<3;j++)
    //for(int j = 0;j<3;j++)
    //    for(int i=0;i<nBins;i++)
    //        AjiCorrelations[i*3+j]->Draw("same");;
    gPad->SetLogx();
    gPad->SetLogy();
    TLegend * legcorr = plotLegend("right_top",  "", 0.9, 0.5, 0., 0., "Autocorrelation",1);
    legcorr->AddEntry(p1CorrelationGibbs, "Gibbs p_{0}", "l");
    legcorr->AddEntry(p2CorrelationGibbs, "Gibbs p_{1}", "l");
    legcorr->AddEntry(p3CorrelationGibbs, "Gibbs p_{2}", "l");
    //legcorr->AddEntry(AjiCorrelations[0], "Gibbs A_{ji}s", "l");
    legcorr->Draw("same");
}

