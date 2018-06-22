#include "GibbsTestClass.h"
#include <TH1F.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TMath.h>
#include "AliMCMCTemplateFitter.h"
#include "AliMCLogLFitter.h"

//_____________________________________________________________________________
GibbsTestClass::GibbsTestClass()
{
   // Constructor
   fRandom = 0;
   fMCTrue = 0x0;
   MeanPosterior0 = 0x0;
   MeanPosterior1 = 0x0;
   MeanPosterior2 = 0x0;
   WidthPosterior0 = 0x0;
   WidthPosterior1 = 0x0;
   WidthPosterior2 = 0x0;
   fML0 = 0x0;
   fML1 = 0x0;
   fML2 = 0x0;
   fNPar=0;
   fBinDef=0x0;
    fNBins = 100;
    fNPar = 3;
    fMCCounts = new Int_t[fNPar];
    ftruePj = new Double_t[fNPar];
    fMCTrue = new TH1D*[fNPar];
    ftruePj[0] = 1.3;
    fMCCounts[0] = 200;
    ftruePj[1] = 1.4;
    fMCCounts[1] = 200;
    ftruePj[2] = 1.5;
    fMCCounts[2] = 200;
}

//_____________________________________________________________________________
GibbsTestClass::~GibbsTestClass()
{
   // Destructor
   if(fRandom) delete fRandom;
   //if(ftruePj) delete[] ftruePj;
   //if (fMCCounts) delete[] fMCCounts;
   //if(fBinDef) delete fBinDef;
   
   //for(int i=0;i<fNPar;i++)
   //    if(fMCTrue[i]) delete fMCTrue[i];
   //if (fMCTrue) delete[] fMCTrue;
}

void GibbsTestClass::Begin(TTree *tree)
{
    
}

//_____________________________________________________________________________
void GibbsTestClass::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   // TString option = GetOption();

   // Histogram
    
    TF1 * fct0 = new TF1("fct0", "TMath::Exp(-TMath::Abs(x/0.2))", -1., 1.);
    TF1 * fct1 = new TF1("fct1", "TMath::Exp(-TMath::Abs(x/0.5))", -1., 1.);
    TF1 * fct2 = new TF1("fct2", "TMath::Exp(-TMath::Abs((x-0.3)/0.2))", -1., 1.);
    
    fBinDef = new TH1D("BinDef", "", fNBins, -1., 1.);
    
    TH1D * H0 =  (TH1D*) fBinDef->Clone("H0");
    TH1D * H1 =  (TH1D*) fBinDef->Clone("H1");
    TH1D * H2 =  (TH1D*) fBinDef->Clone("H2");
    
    for(int i=1;i<=fNBins;i++)
    {
        H0->SetBinContent(i, fct0->Eval(H0->GetXaxis()->GetBinCenter(i)));
        H1->SetBinContent(i, fct1->Eval(H1->GetXaxis()->GetBinCenter(i)));
        H2->SetBinContent(i, fct2->Eval(H2->GetXaxis()->GetBinCenter(i)));
    }
    
    fMCTrue[0]=H0;
    fMCTrue[1]=H1;
    fMCTrue[2]=H2;
    
    delete fct0;
    delete fct1;
    delete fct2;
    
    
   MeanPosterior0 = new TH1D("MeanPosterior0", "", 100, 0., ftruePj[0]*2.);
   MeanPosterior1 = new TH1D("MeanPosterior1", "", 100, 0., ftruePj[1]*2.);
   MeanPosterior2 = new TH1D("MeanPosterior2", "", 100, 0., ftruePj[2]*2.);
   
   WidthPosterior0 = new TH1D("WidthPosterior0", "", 300, 0., ftruePj[0]);
   WidthPosterior1 = new TH1D("WidthPosterior1", "", 300, 0., ftruePj[1]);
   WidthPosterior2 = new TH1D("WidthPosterior2", "", 300, 0., ftruePj[2]);
    
   fML0 = new TH1D("fML0", "", 100, 0., ftruePj[0]*2.);
   fML1 = new TH1D("fML1", "", 100, 0., ftruePj[1]*2.);
   fML2 = new TH1D("fML2", "", 100, 0., ftruePj[2]*2.);
   
   
   fOutput->Add(MeanPosterior0);
   fOutput->Add(MeanPosterior1);
   fOutput->Add(MeanPosterior2);
   fOutput->Add(WidthPosterior0);
   fOutput->Add(WidthPosterior1);
   fOutput->Add(WidthPosterior2);
   fOutput->Add(fML0);
   fOutput->Add(fML1);
   fOutput->Add(fML2);

   // Random number generator
   fRandom = new TRandom3(0);
   
}

//_____________________________________________________________________________
Bool_t GibbsTestClass::Process(Long64_t)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either GibbsTestClass::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
    
    Int_t nIter = 3000;

    TH1D * Template0 = (TH1D*) fBinDef->Clone("Template0");
    TH1D * Template1 = (TH1D*) fBinDef->Clone("Template1");
    TH1D * Template2 = (TH1D*) fBinDef->Clone("Template2");
    
    Int_t nData = fRandom->Poisson(ftruePj[0]*double(fMCCounts[0]) + ftruePj[1]*double(fMCCounts[1]) + ftruePj[2]*double(fMCCounts[2]));
    TH1D * DataTrue = (TH1D*) fBinDef->Clone("DataTrue");
    DataTrue->Add(fMCTrue[0], double(fMCCounts[0])*ftruePj[0]/fMCTrue[0]->Integral());
    DataTrue->Add(fMCTrue[1], double(fMCCounts[1])*ftruePj[1]/fMCTrue[1]->Integral());
    DataTrue->Add(fMCTrue[2], double(fMCCounts[2])*ftruePj[2]/fMCTrue[2]->Integral());
    
    TH1D * DataHist = (TH1D*) fBinDef->Clone("DataHist");
    DataHist->FillRandom(DataTrue, nData);
    Template0->FillRandom(fMCTrue[0], fMCCounts[0]);
    Template1->FillRandom(fMCTrue[1], fMCCounts[1]);
    Template2->FillRandom(fMCTrue[2], fMCCounts[2]);
    
    TH1D ** MCHists = new TH1D*[3];
    MCHists[0]=Template0;
    MCHists[1]=Template1;
    MCHists[2]=Template2;
    
    AliMCMCTemplateFitter * fitClass = new AliMCMCTemplateFitter(3, MCHists, DataHist);
    //if(MakeAjiAutocorr) fitClass->SetPrepareAjiAutocorrelation();
    //fitClass->SetConstantAjiPrior();
    fitClass->SetBurnInRatio(0.5);
    //fitClass->SetMonotonyPrior(0, 9, 10);
    //fitClass->SetMonotonyPrior(1, 9, 10);
    //fitClass->SetMonotonyPrior(2, 12, 13);
    fitClass->Fit(nIter);
    
    MeanPosterior0->Fill(fitClass->ReturnpjMarginal(0)->GetMean());
    MeanPosterior1->Fill(fitClass->ReturnpjMarginal(1)->GetMean());
    MeanPosterior2->Fill(fitClass->ReturnpjMarginal(2)->GetMean());
    
    WidthPosterior0->Fill(fitClass->ReturnpjMarginal(0)->GetRMS());
    WidthPosterior1->Fill(fitClass->ReturnpjMarginal(1)->GetRMS());
    WidthPosterior2->Fill(fitClass->ReturnpjMarginal(2)->GetRMS());
    
    AliMCLogLFitter * MLFitter = new AliMCLogLFitter(3, (TH1**)MCHists, (TH1*)DataHist);
    MLFitter->IfNoMCParticleIsProbablyA(1);
    MLFitter->Fit();
    double MLpar0 = MLFitter->ReturnParameter(0);
    double MLpar1 = MLFitter->ReturnParameter(1);
    double MLpar2 = MLFitter->ReturnParameter(2);
    fML0->Fill(MLpar0);
    fML1->Fill(MLpar1);
    fML2->Fill(MLpar2);
    delete MLFitter;
    
   
   return kTRUE;
}


 //_____________________________________________________________________________
void GibbsTestClass::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   
   MeanPosterior0 = dynamic_cast<TH1D*>(fOutput->FindObject("MeanPosterior0"));
   MeanPosterior1 = dynamic_cast<TH1D*>(fOutput->FindObject("MeanPosterior1"));
   MeanPosterior2 = dynamic_cast<TH1D*>(fOutput->FindObject("MeanPosterior2"));
   WidthPosterior0 = dynamic_cast<TH1D*>(fOutput->FindObject("WidthPosterior0"));
   WidthPosterior1 = dynamic_cast<TH1D*>(fOutput->FindObject("WidthPosterior1"));
   WidthPosterior2 = dynamic_cast<TH1D*>(fOutput->FindObject("WidthPosterior2"));
   fML0 = dynamic_cast<TH1D*>(fOutput->FindObject("fML0"));
   fML1 = dynamic_cast<TH1D*>(fOutput->FindObject("fML1"));
   fML2 = dynamic_cast<TH1D*>(fOutput->FindObject("fML2"));
   
   MeanPosterior0->SetLineColor(kRed+1);
   MeanPosterior1->SetLineColor(kRed+1);
   MeanPosterior2->SetLineColor(kRed+1);
   
   WidthPosterior0->SetLineColor(kRed+1);
   WidthPosterior1->SetLineColor(kBlue+1);
   WidthPosterior2->SetLineColor(kGreen+1);
   
   TLine * p1Truth = new TLine(ftruePj[0], 0., ftruePj[0], 1.2); p1Truth->SetLineWidth(2); p1Truth->SetLineColor(kBlack);
   TLine * p2Truth = new TLine(ftruePj[1], 0., ftruePj[1], 1.2); p2Truth->SetLineWidth(2); p2Truth->SetLineColor(kBlack);
   TLine * p3Truth = new TLine(ftruePj[2], 0., ftruePj[2], 1.2); p3Truth->SetLineWidth(2); p3Truth->SetLineColor(kBlack);
    
   MeanPosterior0->GetXaxis()->SetTitle("p_{0]");
   
   TCanvas *c1 = new TCanvas("MultiCoreTest", "Mean distributions",800,800);
   c1->Divide(2,2);
   c1->cd(1);
   if (fML0) fML0->Draw();
   if (MeanPosterior0) MeanPosterior0->Draw("same");
   p1Truth->Draw("same");
   TLegend * leg = new TLegend(0.1, 0.7, 0.5, 0.9, "");
   leg->AddEntry(fML0, "ML");
   leg->AddEntry(MeanPosterior0, "Mean Posterior");
   leg->Draw("same");
   c1->cd(2);
   if (fML1) fML1->Draw();
   if (MeanPosterior1) MeanPosterior1->Draw("same");
   p2Truth->Draw("same");
   c1->cd(3);
   if (fML2) fML2->Draw();
   if (MeanPosterior2) MeanPosterior2->Draw("same");
   p3Truth->Draw("same");
   c1->cd(4);
   TPaveText *pt = new TPaveText(.05,.1,.95,.8);
   pt->SetFillColor(kWhite);
   pt->AddText("Name : #mu-#mu_{true} : rms : #sqrt{<x-#mu_{true}>^{2}}");
   pt->AddText(Form("ML(0): %.2f : %.2f : %.2f", fML0->GetMean()-ftruePj[0], fML0->GetRMS(), TMath::Sqrt(TMath::Power(fML0->GetMean()-ftruePj[0],2)+TMath::Power(fML0->GetRMS(),2))));
   pt->AddText(Form("Gibbs(0): %.2f : %.2f : %.2f (%.2f)", MeanPosterior0->GetMean()-ftruePj[0], MeanPosterior0->GetRMS(), TMath::Sqrt(TMath::Power(MeanPosterior0->GetMean()-ftruePj[0],2)+TMath::Power(MeanPosterior0->GetRMS(),2)), WidthPosterior0->GetMean()));
   pt->AddText(Form("ML(1): %.2f : %.2f : %.2f", fML1->GetMean()-ftruePj[1], fML1->GetRMS(), TMath::Sqrt(TMath::Power(fML1->GetMean()-ftruePj[1],2)+TMath::Power(fML1->GetRMS(),2))));
   pt->AddText(Form("Gibbs(1): %.2f : %.2f : %.2f (%.2f)", MeanPosterior1->GetMean()-ftruePj[1], MeanPosterior1->GetRMS(), TMath::Sqrt(TMath::Power(MeanPosterior1->GetMean()-ftruePj[1],2)+TMath::Power(MeanPosterior1->GetRMS(),2)), WidthPosterior1->GetMean()));
   pt->AddText(Form("ML(2): %.2f : %.2f : %.2f", fML2->GetMean()-ftruePj[2], fML2->GetRMS(), TMath::Sqrt(TMath::Power(fML2->GetMean()-ftruePj[2],2)+TMath::Power(fML2->GetRMS(),2))));
   pt->AddText(Form("Gibbs(2): %.2f : %.2f : %.2f (%.2f)", MeanPosterior2->GetMean()-ftruePj[2], MeanPosterior2->GetRMS(), TMath::Sqrt(TMath::Power(MeanPosterior2->GetMean()-ftruePj[2],2)+TMath::Power(MeanPosterior2->GetRMS(),2)), WidthPosterior2->GetMean()));
   pt->Draw();
   c1->Update();
   
   TCanvas *c2 = new TCanvas("MultiCoreWidths", "Width distributions",800,800);
   WidthPosterior0->Draw();
   WidthPosterior1->Draw("same");
   WidthPosterior2->Draw("same");
}
