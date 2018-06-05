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
        if(in->GetBinContent(i) > in->GetBinContent(i+1)) in->SetBinContent(i, in->GetBinContent(i+1)*0.99);
    for(int i=n/2+2;i<=n;i++)
        if(in->GetBinContent(i) > in->GetBinContent(i-1)) in->SetBinContent(i, in->GetBinContent(i-1)*0.99); 
}

void TestFitUsingClass(void)
{
    // Simple Test with 3 bins and two sources
    mystyle();
    TRandom3 * rd = new TRandom3(0);
    gRandom->SetSeed(1123);
    cout << "RNG seed: " << gRandom->GetSeed() << endl;
    // Accept-Reject, using Root functionality. A-R w/ function
    
    // function: x^2 * exp(-x)
    
    int nIter = 10000;
    Bool_t MakeAjiAutocorr = kTRUE;
    int rebinY = 40;
    
    
    Double_t p1True = 1.82;
    Int_t n0 = 1350;
    Double_t p2True = 1.55;
    Int_t n1 = 1110;
    Double_t p3True = 2.71;
    Int_t n2 = 894;
    Double_t p4True = 3.51;
    Int_t n3 = 562;
    
    TFile * inFileTreeenhCharm = new TFile("Templates/outputTreeMCenh060.root");
    TH2D * dcaBeauty;
    dcaBeauty = (TH2D*) inFileTreeenhCharm->Get("BeautyMotherCorrelationHalfRAA");
    dcaBeauty->SetName("dcaBeauty");
    TH2D * dcaCharm;
    dcaCharm = (TH2D*) inFileTreeenhCharm->Get("CorrectedCharm3050IP");
    dcaCharm->SetName("dcaCharm");
    dcaCharm->RebinY(rebinY/10);
    dcaBeauty->RebinY(rebinY/10);
    TH1D * Charm = (TH1D*) dcaCharm->ProjectionY("Charm", 6, 15);
    TH1D * Beauty = (TH1D*) dcaBeauty->ProjectionY("Beauty", 6, 15);
    
    
    TFile * inFileTreeMBlargeCent = new TFile("Templates/outputTreeMCMB060.root");
    TH3D * DCAwithCut3dTreeMBlargeCent = (TH3D*) inFileTreeMBlargeCent->Get("TrackDCApTnoCut"); DCAwithCut3dTreeMBlargeCent->SetName("DCAwithCut3dTreeMBlargeCent");
    double EventsTreesMCMBlargeCent = ((TH1D*)inFileTreeMBlargeCent->Get("EventNr"))->Integral();
    DCAwithCut3dTreeMBlargeCent->GetZaxis()->SetRange(3,3);
    TH2D * DCAwithCutTreeMBConversionlargeCent = (TH2D*) DCAwithCut3dTreeMBlargeCent->Project3D("yx");  DCAwithCutTreeMBConversionlargeCent->SetName("DCAwithCutTreeMBConversionlargeCent"); //DCAwithCutTreeMBConversionlargeCent->RebinY(rebinY);
    DCAwithCut3dTreeMBlargeCent->GetZaxis()->SetRange(4,4);
    TH2D * DCAwithCutTreeMBDalitzlargeCent = (TH2D*) DCAwithCut3dTreeMBlargeCent->Project3D("yx");  DCAwithCutTreeMBDalitzlargeCent->SetName("DCAwithCutTreeMBDalitzlargeCent"); //DCAwithCutTreeMBDalitzlargeCent->RebinY(rebinY);
    DCAwithCutTreeMBConversionlargeCent->RebinY(rebinY);
    DCAwithCutTreeMBDalitzlargeCent->RebinY(rebinY);
    
    
    TH1D * Conversion = (TH1D*) DCAwithCutTreeMBConversionlargeCent->ProjectionY("Conversion", 6, 15);
    TH1D * Dalitz = (TH1D*) DCAwithCutTreeMBDalitzlargeCent->ProjectionY("Dalitz", 6, 15);
    
    Charm->Smooth();
    Beauty->Smooth();
    Conversion->Smooth();
    Dalitz->Smooth();
    MakeTruthMonotonous(Charm);
    MakeTruthMonotonous(Beauty);
    MakeTruthMonotonous(Conversion);
    MakeTruthMonotonous(Dalitz);
    SetHistogram(Charm, kGreen+1);
    SetHistogram(Beauty, kRed+1);
    SetHistogram(Conversion, kGray+1);
    SetHistogram(Dalitz, kBlue+1);
    Charm->Scale(1./Charm->Integral());
    Beauty->Scale(1./Beauty->Integral());
    Conversion->Scale(1./Conversion->Integral());
    Dalitz->Scale(1./Dalitz->Integral());
    Int_t nBins = Charm->GetNbinsX();
    double low= Charm->GetXaxis()->GetBinLowEdge(1);
    double high = Charm->GetXaxis()->GetBinLowEdge(nBins+1);
    cout << "4 * " << nBins << " bins. Leading to " << 4*nBins +4 << " unknown free parameters" << endl;
    
    TH1D * BinDefinition = (TH1D*) Conversion->Clone("BinDefinition");
    BinDefinition->Scale(0.);
    BinDefinition->ResetStats();
    
    Int_t nData = rd->Poisson(p1True*double(n0) + p2True*double(n1) + p3True*double(n2) + p4True*double(n3));
    TH1D * DataTrue = (TH1D*) BinDefinition->Clone("DataTrue");
    DataTrue->Add(Charm, double(n0)*p1True/Charm->Integral());
    DataTrue->Add(Beauty, double(n1)*p2True/Beauty->Integral());
    DataTrue->Add(Conversion, double(n2)*p3True/Conversion->Integral());
    DataTrue->Add(Dalitz, double(n3)*p4True/Dalitz->Integral());
    
    TH1D * DataHist = (TH1D*) BinDefinition->Clone("DataHist");
    DataHist->FillRandom(DataTrue, nData);
    TH1D * Template0 = (TH1D*) BinDefinition->Clone("Template0");
    Template0->FillRandom(Charm, n0);
    TH1D * Template1 = (TH1D*) BinDefinition->Clone("Template1");
    Template1->FillRandom(Beauty, n1);
    TH1D * Template2 = (TH1D*) BinDefinition->Clone("Template2");
    Template2->FillRandom(Conversion, n2);
    TH1D * Template3 = (TH1D*) BinDefinition->Clone("Template3");
    Template3->FillRandom(Dalitz, n3);
    cout << "Maximum template1: " << Template0->GetMaximum() << " template2: " << Template1->GetMaximum() << " template3: " << Template2->GetMaximum() << " template4: " << Template3->GetMaximum() << " data: " << DataHist->GetMaximum() <<  endl;
    cout << "Largest bin : " << Charm->GetMaximumBin() << endl;
    
    
    // Now try with class
    TH1D ** MCHists = new TH1D*[4];
    MCHists[0]=Template0;
    MCHists[1]=Template1;
    MCHists[2]=Template2;
    MCHists[3]=Template3;
    AliMCMCTemplateFitter * fitClass = new AliMCMCTemplateFitter(4, MCHists, DataHist);
    if(MakeAjiAutocorr) fitClass->SetPrepareAjiAutocorrelation();
    fitClass->SetBurnInRatio(0.2);
    //fitClass->SetMonotonyPrior(0, 48, 50);
    //fitClass->SetMonotonyPrior(1, 48, 50);
    //fitClass->SetMonotonyPrior(2, 46, 48);
    //fitClass->SetMonotonyPrior(3, 48, 50);
    fitClass->Fit(nIter);
    TH1D * p1GibbsClass = fitClass->ReturnpjMarginal(0);
    TH1D * p2GibbsClass = fitClass->ReturnpjMarginal(1);
    TH1D * p3GibbsClass = fitClass->ReturnpjMarginal(2);
    TH1D * p4GibbsClass = fitClass->ReturnpjMarginal(3);
    TH2D * pPosteriorGibbs = fitClass->ReturnpjMarginal(0,1);
    //TH2D * pPosteriorGibbs12 = fitClass->ReturnpjMarginal(1,2);
    TH2D * pPosteriorGibbs23 = fitClass->ReturnpjMarginal(2,3);
    //pPosteriorGibbs->SetAxisRange(0.,2.,"X");
    //pPosteriorGibbs->SetAxisRange(0.,2.,"Y");
    p1GibbsClass->RebinX(5);
    p2GibbsClass->RebinX(5);
    p3GibbsClass->RebinX(5);
    p4GibbsClass->RebinX(5);
    
    
    cout << "Calculating correlations." << endl;
    
    int MaxCorr = 1000;
    
    TH1D * p1CorrelationGibbs = fitClass->GetpjAutocorrelation(0, MaxCorr);
    TH1D * p2CorrelationGibbs = fitClass->GetpjAutocorrelation(1, MaxCorr);
    TH1D * p3CorrelationGibbs = fitClass->GetpjAutocorrelation(2, MaxCorr);
    TH1D * p4CorrelationGibbs = fitClass->GetpjAutocorrelation(3, MaxCorr);
    
    TH1D ** AjiCorrelations = new TH1D*[4*nBins];
    
    TH2D * A0iMarginal = fitClass->ReturnAjiMarginals(0);
    TH2D * A1iMarginal = fitClass->ReturnAjiMarginals(1);
    TH2D * A2iMarginal = fitClass->ReturnAjiMarginals(2);
    TH2D * A3iMarginal = fitClass->ReturnAjiMarginals(3);
    if(MakeAjiAutocorr)
    for(int j = 0;j<4;j++)
        for(int i=0;i<nBins;i++)
        {
            AjiCorrelations[i*4+j] = fitClass->GetAjiAutocorrelation(j, i, MaxCorr);
            AjiCorrelations[i*4+j]->SetLineColor(kBlack);
        }
    
    
    cout << "Done." << endl;
    
    cout << "Now checking with the ML fit." << endl;
    AliMCLogLFitter * MLFitter = new AliMCLogLFitter(4, (TH1**)MCHists, (TH1*)DataHist);
    //Fitter->SetUseChiSq();
    MLFitter->IfNoMCParticleIsProbablyA(1);
    //Fitter->SetFitRange(-0.10, 0.10);
    MLFitter->Fit();
    double MLpar0 = MLFitter->ReturnParameter(0);
    double MLpar1 = MLFitter->ReturnParameter(1);
    double MLpar2 = MLFitter->ReturnParameter(2);
    double MLpar3 = MLFitter->ReturnParameter(3);
    TMarker * MLMarker01 = new TMarker(MLpar0, MLpar1, 20); MLMarker01->SetMarkerColor(kGreen+2);
    TMarker * MLMarker23 = new TMarker(MLpar2, MLpar3, 20); MLMarker23->SetMarkerColor(kGreen+2);
    cout << "end of ML fit" << endl;
    
    TH1D * BlankHisto = new TH1D("bg", "", 100, low, high);
    BlankHisto->SetLineColor(kWhite);
    BlankHisto->SetAxisRange(0.1, DataHist->GetMaximum()*1.2,"Y");
    BlankHisto->GetXaxis()->SetTitle("x");
    BlankHisto->GetYaxis()->SetTitle("counts");
    SetHistogram(DataHist, kBlack);
    SetHistogram(Template0, kGreen+1);
    SetHistogram(Template1, kRed+1);
    SetHistogram(Template2, kGray+1);
    SetHistogram(Template3, kBlue+1);
    TH1D * CharmInputCopy = (TH1D*) Template0->Clone("CharmInputCopy");
    TH1D * BeautyInputCopy = (TH1D*) Template1->Clone("BeautyInputCopy");
    TH1D * ConversionInputCopy = (TH1D*) Template2->Clone("ConversionInputCopy");
    TH1D * DalitzInputCopy = (TH1D*) Template3->Clone("DalitzInputCopy");
    Template0->Scale(p1GibbsClass->GetMean());
    Template1->Scale(p2GibbsClass->GetMean());
    Template2->Scale(p3GibbsClass->GetMean());
    Template3->Scale(p4GibbsClass->GetMean());
    TH1D * Sum = (TH1D*) BinDefinition->Clone("Sum");
    Sum->Add(Template0);
    Sum->Add(Template1);
    Sum->Add(Template2);
    Sum->Add(Template3);
    Sum->SetLineWidth(2);Sum->SetLineColor(kOrange+1);
    p1GibbsClass->Scale(1./p1GibbsClass->Integral("width"));
    p2GibbsClass->Scale(1./p2GibbsClass->Integral("width"));
    p3GibbsClass->Scale(1./p3GibbsClass->Integral("width"));
    p4GibbsClass->Scale(1./p4GibbsClass->Integral("width"));
    SetHistogram(p1GibbsClass, kGreen+1);
    SetHistogram(p2GibbsClass, kRed+1);
    SetHistogram(p3GibbsClass, kGray+1);
    SetHistogram(p4GibbsClass, kBlue+1);
    SetHistogram(p1CorrelationGibbs, kGreen+1);
    SetHistogram(p2CorrelationGibbs, kRed+1);
    SetHistogram(p3CorrelationGibbs, kGray+1);
    SetHistogram(p4CorrelationGibbs, kBlue+1);
    p2CorrelationGibbs->SetLineStyle(2);
    p1CorrelationGibbs->SetAxisRange(0.07, 1.1, "Y");
    
    TMarker * Truth = new TMarker(p1True, p2True, 20);
    TMarker * Truth12 = new TMarker(p2True, p3True, 20);
    TMarker * Truth23 = new TMarker(p3True, p4True, 20);
    
    TH2D * cont12 = fitClass->GetCredibleArea(1, 2, 0.95);
    TH2D * cont23 = fitClass->GetCredibleArea(2, 3, 0.95);
    TH2D * cont01 = fitClass->GetCredibleArea(0, 1, 0.95);
    
    double MarginalMaximum = TMath::Max(TMath::Max(p1GibbsClass->GetMaximum(), p2GibbsClass->GetMaximum()), p3GibbsClass->GetMaximum());
    p2GibbsClass->SetAxisRange(0., MarginalMaximum*1.2, "Y");
    
    TLine * p1Truth = new TLine(p1True, 0., p1True, MarginalMaximum*0.1); p1Truth->SetLineWidth(2); p1Truth->SetLineColor(kGreen+1);
    TLine * p2Truth = new TLine(p2True, 0., p2True, MarginalMaximum*0.1); p2Truth->SetLineWidth(2); p2Truth->SetLineColor(kRed+1);
    TLine * p3Truth = new TLine(p3True, 0., p3True, MarginalMaximum*0.1); p3Truth->SetLineWidth(2); p3Truth->SetLineColor(kGray+1);
    TLine * p4Truth = new TLine(p4True, 0., p4True, MarginalMaximum*0.1); p4Truth->SetLineWidth(2); p4Truth->SetLineColor(kBlue+1);
    
    cout << "Uncertainties:\nCharm: " << p1GibbsClass->GetRMS()
        << "\nBeauty: " << p2GibbsClass->GetRMS()
        << "\nConversion: " << p3GibbsClass->GetRMS()
        << "\nDalitz: " << p4GibbsClass->GetRMS() << endl;
        
    
    
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
    Template3->Draw("sameehist");
    gPad->SetLogy();
    TLegend * leg = plotLegend("right_top",  "", 0.9, 0.7, -0.05, 0., "",1);
    leg->AddEntry(DataHist, "Data", "l");
    leg->AddEntry(Template0, "Template0*p_{0}", "l");
    leg->AddEntry(Template1, "Template1*p_{1}", "l");
    leg->AddEntry(Template2, "Template2*p_{2}", "l");
    leg->AddEntry(Template3, "Template3*p_{3}", "l");
    leg->AddEntry(Sum, "Sum of scaled dist", "l");
    leg->Draw("same");
    cv->cd(2);
    pPosteriorGibbs->Draw("colz");
    SetSquare2DHist((TH2*)pPosteriorGibbs, (TPad*)gPad);
    cont01->Draw("samecont2");
    Truth->Draw("same");
    MLMarker23->Draw("same");
    TLegend * leg01 = plotLegend("right_top",  "", 1.2, 0.7, -0.2, 0., Form("Gibbs, %d iterations", nIter),1);
    leg01->AddEntry(Truth, "True value", "p");
    leg01->AddEntry(cont01, "95% credible area", "l");
    leg01->AddEntry(MLMarker23, "ML", "p");
    leg01->Draw("same");
    cv->cd(3);
    pPosteriorGibbs23->Draw("colz");
    SetSquare2DHist((TH2*)pPosteriorGibbs23, (TPad*)gPad);
    Truth23->Draw("same");
    MLMarker23->Draw("same");
    cont23->Draw("samecont2");
    cv->cd(4);
    Charm->Draw("hist");
    SetRect1DHist((TH1*)Charm, (TPad*)gPad);
    Beauty->Draw("SAMEhist");
    Conversion->Draw("samehist");
    Dalitz->Draw("samehist");
    gPad->SetLogy();
    cv->cd(5);
    p2GibbsClass->SetAxisRange(0., TMath::Max(TMath::Max(p1GibbsClass->GetMaximum(), p2GibbsClass->GetMaximum()), p3GibbsClass->GetMaximum())*1.2, "Y");
    p2GibbsClass->Draw("hist");
    p1GibbsClass->Draw("samehist");
    p3GibbsClass->Draw("samehist");
    p4GibbsClass->Draw("samehist");
    SetRect1DHist((TH1*)p2GibbsClass, (TPad*)gPad);
    p2Truth->Draw("same");
    p1Truth->Draw("same");
    p3Truth->Draw("same");
    p4Truth->Draw("same");
    TLegend * legp2 = plotLegend("right_top",  "", 1.1, 0.8, -0.1, 0., "marginals with true values",1);
    legp2->AddEntry(p1GibbsClass, "p_{0}", "l");
    legp2->AddEntry(p2GibbsClass, "p_{1}", "l");
    legp2->AddEntry(p3GibbsClass, "p_{2}", "l");
    legp2->AddEntry(p4GibbsClass, "p_{3}", "l");
    legp2->Draw("same");
    cv->cd(6);
    p1CorrelationGibbs->Draw("");
    SetRect1DHist((TH1*)p1CorrelationGibbs, (TPad*)gPad);
    p2CorrelationGibbs->Draw("same");
    p3CorrelationGibbs->Draw("same");
    p4CorrelationGibbs->Draw("same");
    if(MakeAjiAutocorr)
    for(int j = 0;j<4;j++)
        for(int i=0;i<nBins;i++)
            AjiCorrelations[i*4+j]->Draw("same");;
    gPad->SetLogx();
    gPad->SetLogy();
    TLegend * legcorr = plotLegend("right_top",  "", 0.9, 0.6, 0., 0., "Autocorrelation",1);
    legcorr->AddEntry(p1CorrelationGibbs, "p_{0}", "l");
    legcorr->AddEntry(p2CorrelationGibbs, "p_{1}", "l");
    legcorr->AddEntry(p3CorrelationGibbs, "p_{2}", "l");
    legcorr->AddEntry(p4CorrelationGibbs, "p_{3}", "l");
    //legcorr->AddEntry(AjiCorrelations[0], "Gibbs A_{ji}s", "l");
    legcorr->Draw("same");
    
    
    TCanvas * cv2 = new TCanvas("BeautyAjiMarginals", "", 1400, 1000);
    cv2->Divide(2,2);
    cv2->cd(1);
    A0iMarginal->Draw("colz");
    CharmInputCopy->Draw("histesame");
    gPad->SetLogy();
    gPad->SetLogz();
    cv2->cd(2);
    A1iMarginal->Draw("colz");
    BeautyInputCopy->Draw("histesame");
    gPad->SetLogy();
    gPad->SetLogz();
    cv2->cd(3);
    A2iMarginal->Draw("colz");
    ConversionInputCopy->Draw("histesame");
    gPad->SetLogy();
    gPad->SetLogz();
    cv2->cd(4);
    A3iMarginal->Draw("colz");
    DalitzInputCopy->Draw("histesame");
    gPad->SetLogy();
    gPad->SetLogz();
    
}

