{
    gSystem->Load("libMathMore.so");
    gROOT->LoadMacro("AliMCMCTemplateFitter.cxx+");
    gROOT->LoadMacro("AliMCLogLFitter.cxx+");
    //gROOT->LoadMacro("LargerBinNumber.C");
    //gROOT->LoadMacro("ThreeSources.C");
    //gROOT->LoadMacro("TestUsingIPDistributions.C");
    gROOT->LoadMacro("TestUsingExponentialTemplates.C");
    //gROOT->LoadMacro("TestFitUsingClass2.C");
    TestFitUsingClass();
}
