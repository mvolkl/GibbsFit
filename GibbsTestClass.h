#ifndef GibbsTestClass_h
#define GibbsTestClass_h

#include <TSelector.h>

class TH1F;
class TH1D;
class TRandom;
class GibbsTestClass : public TSelector {
public :

   // Define members here
   TH1D ** fMCTrue;
   Int_t fNPar;
   Int_t fNBins;
   TH1D * fBinDef;
   Double_t * ftruePj;
   Int_t * fMCCounts;
   
   // Output Objects
   TH1D  * MeanPosterior0; //
   TH1D  * MeanPosterior1;
   TH1D  * MeanPosterior2;
   TH1D  * WidthPosterior0; //
   TH1D  * WidthPosterior1;
   TH1D  * WidthPosterior2;
   TH1D  * fML0; //
   TH1D  * fML1;
   TH1D  * fML2;
   TRandom *fRandom;  //! Random number generator

   GibbsTestClass();
   virtual ~GibbsTestClass();
   virtual Int_t   Version() const { return 1; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual Bool_t  Process(Long64_t entry);
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate(){}
   virtual void    Terminate();

   ClassDef(GibbsTestClass,1);
};
#endif
