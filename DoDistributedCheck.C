{
    gSystem->Load("libMathMore.so");
    gROOT->LoadMacro("AliMCMCTemplateFitter.cxx+");
    gROOT->LoadMacro("AliMCLogLFitter.cxx+");
    //gROOT->LoadMacro("TestUsingExponentialTemplates.C");
    TProof *plite= TProof::Open("lite://", "workers=4");
    plite->Process("GibbsTestClass.C+", 1000);
}
