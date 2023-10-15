// main.cc --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 一 9月 21 20:06:21 2020 (+0800)
// Last-Updated: 一 1月 16 11:43:06 2023 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 6
// URL: http://wuhongyi.cn 

#include "DataAnalysis.hh"
#include "UserDefine.hh"

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TRint.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const char gVERSION[] = "Version: DataAnalysis - 2023.01.16";

int main(int argc, char* argv[])
{
  std::cout<<gVERSION<<std::endl;

  //TRint *theApp = new TRint("Rint",&argc,argv);
    
  if(argc < 2)
    {
      std::cout << "eg. " << argv[0] << "  [run number]" <<std::endl;
      return 1;
    }

  TString RunNumber(argv[1]);
  int run = RunNumber.Atoi();
   
  Long64_t startentry = -1, stopentry = -1;
  if (argc >= 4){
    TString startNumber(argv[2]);
    startentry = startNumber.Atoi();

    TString stopNumber(argv[3]);
    stopentry = stopNumber.Atoi();
  }

  TChain *ipt = new TChain("tree");
#ifdef WAVEFORM
  ipt->Add(TString::Format("%s%s_%04d_W%d_wave.root", RAWFILEPATH, RAWFILENAME, run, EVENTTIMEWINDOWSWIDTH).Data());
  ipt->Add(TString::Format("%s%s_%04d_W%d_wave_*.root", RAWFILEPATH, RAWFILENAME, run, EVENTTIMEWINDOWSWIDTH).Data());
#else
  ipt->Add(TString::Format("%s%s_%04d_W%d_notwave.root",RAWFILEPATH, RAWFILENAME, run, EVENTTIMEWINDOWSWIDTH).Data());
#endif  
    
  TFile *ipf_wave = new TFile(Form("%s%s_%04d_W%d.root", WAVEFILEPATH, WAVEFILENAME, run, EVENTTIMEWINDOWSWIDTH));

  TString opfname = TString::Format("%s%s_%04d_W%d.root", ROOTFILEPATH, ROOTFILENAME, run, EVENTTIMEWINDOWSWIDTH);
  if (argc >= 4){
      if ( startentry > ipt->GetEntries() && stopentry > ipt->GetEntries() )
          return 0;
      
      opfname = TString::Format("%s%s_%04d_W%d_%lld_%lld.root", ROOTFILEPATH, ROOTFILENAME, run, EVENTTIMEWINDOWSWIDTH, startentry, stopentry);
  }
    
  auto opf = new TFile(opfname.Data(), "RECREATE");
  auto opt = new TTree("tree", "tree");

  auto it = new DataAnalysis(ipt, ipf_wave);
  it->Loop(opf, opt, startentry, stopentry);

  opt->Write();
  opf->Close();
  ipf_wave->Close();
    
  //theApp->Run();
  //delete theApp;
  
  return 0;
}

// 
// main.cc ends here
