#include "mystyle.h"

#include <TPaveText.h>
#include <TColor.h>

void SetSquare2DHist(TH2* hist, TPad * pad)
{
  pad->SetRightMargin(0.14);
  pad->SetBottomMargin(0.15);
  pad->SetLeftMargin(0.11);
  pad->SetTopMargin(0.1);
  hist->GetXaxis()->SetTitleSize(0.04);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetTitleOffset(1.0);
  hist->GetYaxis()->SetTitleSize(0.04);
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetZaxis()->SetTitleSize(0.04);
  hist->GetZaxis()->SetLabelSize(0.04);
  hist->GetZaxis()->SetTitleOffset(1.0);
  pad->SetFrameLineWidth(2);
}

void SetRect1DHist(TH1* hist, TPad * pad)
{
  hist->SetTitle("");
  pad->SetRightMargin(0.05);
  pad->SetBottomMargin(0.15);
  pad->SetLeftMargin(0.10);
  pad->SetTopMargin(0.1);
  pad->SetFrameLineWidth(2);
  
  hist->GetXaxis()->SetTitleSize(0.04);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetTitleSize(0.04);
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetZaxis()->SetTitleSize(0.04);
  hist->GetZaxis()->SetLabelSize(0.04);
  hist->GetZaxis()->SetTitleOffset(1.0);  
  hist->SetLineWidth(2);
}

void SetRect2DHist(TH1* hist, TPad * pad)
{
  hist->SetTitle("");
  pad->SetRightMargin(0.1);
  pad->SetBottomMargin(0.15);
  pad->SetLeftMargin(0.10);
  pad->SetTopMargin(0.1);
  pad->SetFrameLineWidth(2);
  
  hist->GetXaxis()->SetTitleSize(0.04);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetTitleSize(0.04);
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetZaxis()->SetTitleSize(0.04);
  hist->GetZaxis()->SetLabelSize(0.04);
  hist->GetZaxis()->SetTitleOffset(1.0);  
  hist->SetLineWidth(2);
}


void WriteThisThesis(double x, double y, double factor=1.0)
{
   TPaveText *pt = new TPaveText(x,y,x+.2*factor,y+0.1*factor, "ndc");
   pt->AddText("This Thesis");
   pt->SetFillStyle(0);
   pt->SetBorderSize(0);
   //pt->SetTextSize(.1);
   //pt->AddText("They are added to the pave using the AddText method.");
   //pt->AddLine(.0,.5,1.,.5);
   //pt->AddText("Even complex TLatex formulas can be added:");
   //pt->AddText("F(t) = #sum_{i=-#infty}^{#infty}A(i)cos#[]{#frac{i}{t+i}}");
   pt->Draw();
}

void SetBox(TBox * box, int color)
{
  box->SetFillStyle(0);
  box->SetLineWidth(2);
  box->SetLineColor(color);
  
}

TPad **  MakeRatioPlot(TH1D * UpperFunction, TH1D * LowerFunction, TPad * pad)
{
  TPad ** pads = new TPad*[2];
  
  double factorUpper = 4./3.;
  double factorLower = 4.;
  double LabelSize = 0.04;  // 0.05
  double TitleSize = 0.04; // 0.06
  double topOffset = 0.05;
  double tickLength = 0.03;
  
  TPad * lowerPad = new TPad ("Ratio", "", 0.0, 0.00, 1., 0.255);
  lowerPad->SetFrameLineWidth(2);
  lowerPad->Draw();
  lowerPad->cd();
  //setOPT_canvas((TCanvas*)lowerPad);
  lowerPad->SetFillStyle(4000);
  lowerPad->SetTopMargin(0.01);
  lowerPad->SetBottomMargin(0.4);
  LowerFunction->SetTitle("");
  //setOPT_hists((TH1F*) LowerFunction,"d_{0} #times sgn(charge #times field) (cm)","Data/Fit ",510,20,1.0, kBlack);
  LowerFunction->GetYaxis()->SetNdivisions(3);
  LowerFunction->GetYaxis()->SetTitleSize(TitleSize*factorLower*1.2);
  LowerFunction->GetYaxis()->SetTitleOffset(0.05*factorLower);
  LowerFunction->GetYaxis()->SetTickLength(tickLength);
  //LowerFunction->GetYaxis()->SetNdivisions(3);
  LowerFunction->GetXaxis()->SetTitleSize(TitleSize*factorLower*1.2);
  LowerFunction->GetXaxis()->SetTitleOffset(0.3*factorLower/1.4);
  LowerFunction->GetYaxis()->SetLabelSize(LabelSize*factorLower);
  LowerFunction->GetXaxis()->SetLabelSize(LabelSize*factorLower);
  LowerFunction->GetXaxis()->SetTickLength(tickLength*factorLower);
  //lowerPad->SetFrameLineWidth(5);
  //LowerFunction->Draw("hist");
  LowerFunction->Draw("hist");
  TF1 * one = new TF1("one", "1", LowerFunction->GetXaxis()->GetBinLowEdge(1), LowerFunction->GetXaxis()->GetBinUpEdge(LowerFunction->GetNbinsX())); one->SetLineColor(kBlack); one->SetLineWidth(1); one->SetLineStyle(2);
  one->Draw("SAME");
  pad->cd(0);
  
  
  
  TPad * upperPad = new TPad ("Distributions", "", 0.0, 0.25, 1., 1.);
  upperPad->SetFrameLineWidth(2);
  upperPad->Draw();
  upperPad->cd();
  upperPad->SetBottomMargin(0.0);
  upperPad->SetTopMargin(topOffset* factorUpper);
  //setOPT_canvas((TCanvas*)upperPad);
  upperPad->SetFillStyle(4000);
  
  UpperFunction->GetYaxis()->SetTitleSize(TitleSize*factorUpper*1.2);
  UpperFunction->GetYaxis()->SetTitleOffset(0.5*factorUpper);
  UpperFunction->GetYaxis()->SetTickLength(tickLength);
  UpperFunction->GetXaxis()->SetTitleSize(TitleSize*factorUpper*1.2);
  UpperFunction->GetXaxis()->SetTitleOffset(0.3*factorUpper);
  UpperFunction->GetYaxis()->SetLabelSize(LabelSize*factorUpper);
  UpperFunction->GetXaxis()->SetLabelSize(LabelSize*factorUpper);
  UpperFunction->GetXaxis()->SetTickLength(tickLength*factorUpper);
  //UpperFunction->SetAxisRange(2, UpperFunction->GetMaximum()*4., "Y");
  upperPad->SetLogy();
  UpperFunction->Draw("ep");
  
  /*pad->cd(0);
  TLine * line = new TLine(0.01, 0.25, 0.9, 0.25);
  line->SetLineWidth(20.);
  line->Draw("SAME");
  line->Pop();*/
  
  
  lowerPad->Pop();
  
  
  pads[0]=upperPad;
  pads[1]=lowerPad;
  return pads;
}

void SetBetterRainbowPalette(void)
{
  gStyle->SetPalette(55);
}

void SetLinearPalette(void) // from negative to positive values
{
  double stops[5] = {0., 0.4, 0.5, 0.6 ,1.0};
  double r[5] = {0., 0.5,  1., 0.8, 0.5};
  double g[5] = {0., 0.5,  1., 0.5, 0.0};
  double b[5] = {0.5, 0.8,1., 0.5, 0.};
  TColor::CreateGradientColorTable(5, &stops[0], &r[0],&g[0],&b[0],1000);
  gStyle->SetNumberContours(100);
}

void SetScalePalette(void) // from 0 to positive values
{
  double stops[5] = {0., 0.1, 0.33, 0.66, 1.0};
  double r[5] = {1., 0.8, 0.5, 0., 0.0};
  double g[5] = {1., 0.8, 1., 0., 0.0};
  double b[5] = {1.0, 0.8, 0.5, 1.0, 0.};
  TColor::CreateGradientColorTable(5, &stops[0], &r[0],&g[0],&b[0],1000);
  gStyle->SetNumberContours(100);
}

void SetScalePalette2(void) // from 0 to positive values
{
  double stops[3] = {0., 0.5, 1.0};
  double r[3] = {1., 0., 0.0};
  double g[3] = {1., 0., 0.0};
  double b[3] = {1.0, 1.0, 0.};
  TColor::CreateGradientColorTable(3, &stops[0], &r[0],&g[0],&b[0],1000);
  gStyle->SetNumberContours(100);
}

void SetZAxisLinear(TH2D* in)
{
  double max = TMath::Max(in->GetMaximum(), -in->GetMinimum());
  in->SetAxisRange(-max, max*1.01, "Z");
}

void UseGrid(TH1 * in)
{
  gPad->SetGridx();
  gPad->SetGridy();
}
