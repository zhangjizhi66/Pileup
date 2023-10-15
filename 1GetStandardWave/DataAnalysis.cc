// event.cc --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 一 9月 21 16:28:37 2020 (+0800)
// Last-Updated: 四 4月 27 22:57:22 2023 (+0800)
//	     By: Hongyi Wu(吴鸿毅)
//     Update #: 123
// URL: http://wuhongyi.cn

#include "DataAnalysis.hh"

#include "UserDefine.hh"
#include "TMath.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TDirectoryFile.h"
#include "TROOT.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DataAnalysis::Init()
{
    if (ipt == NULL) return;
    ipt->SetBranchAddress("event", &br_event);
    
    for (int idet = 0; idet < ndet; idet++)
        for (int id = 0; id < nid; id++){
            if ( !valid(idet,id) ) continue;
            
            h2fwhm[idet][id] = (TH2D*)ipf_check->Get(Form("h2fwhm_det%d_id%d",idet,id));
            h2fwhm[idet][id]->ProjectionX(Form("tpjx_det%d_id%d",idet,id));
            TH1F *tpjx = (TH1F*)gROOT->FindObject(Form("tpjx_det%d_id%d",idet,id));
            TH1F *hgrx = (TH1F*)tpjx->Clone("hgrx");
            hgrx->Reset();
            
            gROOT->Macro(Form("%scut_%04d_%d_%d.C",CHECKFILEPATH,run,idet,id));
            cutg[idet][id] = (TCutG*)gROOT->FindObject(Form("cut_%04d_%d_%d",run,idet,id));
            
            grfwhm_raw[idet][id] = (TGraph*)ipf_check->Get(Form("grfwhm_det%d_id%d",idet,id));
            grfwhm[idet][id] = new TGraph;
            grfwhm[idet][id]->SetName(Form("grfwhm_det%d_id%d",idet,id));
            for (int ipnt = 0; ipnt < grfwhm_raw[idet][id]->GetN(); ipnt++)
                if ( cutg[idet][id]->IsInside(grfwhm_raw[idet][id]->GetPointX(ipnt), grfwhm_raw[idet][id]->GetPointY(ipnt)) ){
                    double x = grfwhm_raw[idet][id]->GetPointX(ipnt);
                    double y = grfwhm_raw[idet][id]->GetPointY(ipnt);
                    int xbin = hgrx->FindBin(x);
                    int xbin0 = hgrx->FindBin(refbin);
                    
                    if ( hgrx->GetBinContent(xbin) >= tpjx->GetBinContent(xbin0) ) continue;
                    hgrx->Fill(x);
                    
                    grfwhm[idet][id]->SetPoint(grfwhm[idet][id]->GetN(), x, y);
                }
            
            fpol2[idet][id] = new TF1(Form("fpol2_det%d_id%d",idet,id), "pol2", 0, maxoverflow);
            grfwhm[idet][id]->Fit(fpol2[idet][id], "Q ROB");
            grfwhm[idet][id]->Write();
        }
  
    // init wave
    for (int idet = 0; idet < ndet; idet++)
        for (int id = 0; id < nid; id++){
            if ( !valid(idet,id) ) continue;
            
            standardwave[idet][id] = new TGraph;
            standardwave[idet][id]->SetName(Form("standardwave_det%d_id%d",idet,id));
            standardwave[idet][id]->SetUniqueID(0);

            for (int ipnt = 0; ipnt < 2*length; ipnt++)
                standardwave[idet][id]->SetPoint(ipnt, ipnt, 0);
        }
}

void DataAnalysis::Loop(TTree *opt_)
{
    if (opt_ == NULL) return;

    opt = opt_;
    BranchOpt();
    
    clock_t start = clock(), stop = clock();
  
    Long64_t startentry = 0, stopentry = ipt->GetEntries();
    for (Long64_t jentry = startentry; jentry < stopentry; jentry++){
        ipt->GetEntry(jentry);
       
        // vref.clear();
        vdet.clear();

        for (int i = 0; i < int(br_event->size()); i++){
            int flag = (*br_event)[i].det;
            det.det = (*br_event)[i].det;
            det.id = (*br_event)[i].id;
            det.raw = (*br_event)[i].raw;
            det.e = (*br_event)[i].e;
            det.ts = (*br_event)[i].ts;
            det.pileup = (*br_event)[i].pileup;
            det.outofr = (*br_event)[i].outofr;

            det.subts = 0.0;

            if ( (*br_event)[i].sr == 250 )
                det.ts = 8*det.ts;
            else if ( (*br_event)[i].sr == 100 || (*br_event)[i].sr == 500 )
                det.ts = 10*det.ts;

            if ( (*br_event)[i].sr == 500 )
                det.subts = ((*br_event)[i].cfds-1+(*br_event)[i].cfd/8192.0)*2.0;
            else if ( (*br_event)[i].sr == 250 )
                det.subts = ((*br_event)[i].cfd/16384.0-(*br_event)[i].cfds)*4.0;
            else if ( (*br_event)[i].sr == 100 )
                det.subts = ((*br_event)[i].cfd/32768.0)*10.0;

#ifdef WAVEFORM
            det.wave.clear();
            if ( (*br_event)[i].ltra > 0 ){     // !!!!!!
                det.wave.assign((*br_event)[i].data.begin(),(*br_event)[i].data.end());
                // det.wave.size()    det.wave[i]
            }
#endif
            //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

            if ( flag == 0 ){     // ref for beam bursts
                vref.clear();
                vref.push_back(det);
            }
            else
                vdet.push_back(det);
        }
      
        for (int ihit = 0; ihit < int(vdet.size()); ihit++){
            if ( int(vdet[ihit].wave.size()) <= 0 ) continue;
            
            if ( !valid(vdet[ihit].det, vdet[ihit].id) ) continue;
            
            // get baseline
            double base = 0;
            for (int ipnt = 0; ipnt < nbase; ipnt++)
                base += (double)vdet[ihit].wave[ipnt] / nbase;
            
            // fastfilter
            std::vector<int> vtrig;
            int last1fastfilter = 0, last2fastfilter = 0;
            int thres = 0;
            for (int ipnt = 0; ipnt < int(vdet[ihit].wave.size()); ipnt++){
                int fastfilter = 0;
                if ( ipnt >= L && ipnt <= int(vdet[ihit].wave.size())-L ){
                    for (int jpnt = 0; jpnt < L; jpnt++)
                        fastfilter += vdet[ihit].wave[ipnt+jpnt] - vdet[ihit].wave[ipnt-L+jpnt];
                    fastfilter /= L;
                }

                if ( ipnt < nbase ){
                    if ( fastfilter > thres )
                        thres = fastfilter;
                }
                else if ( ipnt == nbase )
                    thres *= 2;
                else {
                    if ( fastfilter >= thres && last1fastfilter < thres && last2fastfilter < thres )
                        vtrig.push_back(ipnt);
                }
                
                last2fastfilter = last1fastfilter;
                last1fastfilter = fastfilter;
            }

            // get standard wave
            for (int itrig = 0; itrig < int(vtrig.size()); itrig++){
                
                // pileup cut
                if ( itrig > 0 && vtrig[itrig] - vtrig[itrig-1] < length ) continue;
                if ( itrig < int(vtrig.size())-1 && vtrig[itrig+1] - vtrig[itrig] < length ) continue;
                if ( vtrig[itrig] < pileupzone ) continue;
                
                // peaks
                int xmax = -1;
                double ymax = -1;
                for (int ipnt = vtrig[itrig] + trise/2; ipnt < int(vdet[ihit].wave.size()); ipnt++){
                    if ( vdet[ihit].wave[ipnt] - base > ymax ){
                        ymax = vdet[ihit].wave[ipnt] - base;
                        xmax = ipnt;
                    }
                    else if (ipnt - xmax > 5)
                        break;
                }

                // position cut
                if ( xmax < length || xmax > int(vdet[ihit].wave.size()) - length ) continue;
                
                // height cut
                if ( ymax < minheight || ymax > maxheight ) continue;
                
                // fill standard wave
                double wavefwhm = fpol2[vdet[ihit].det][vdet[ihit].id]->Eval(ymax);
                TGraph *waveshrink = new TGraph;
                for (int ipnt = 0; ipnt < 2*length; ipnt++)
                    waveshrink->SetPoint(ipnt, length+(ipnt-length)*standardfwhm/wavefwhm, (vdet[ihit].wave[xmax-length+ipnt] - base) / ymax);
                for (int ipnt = 0; ipnt < 2*length; ipnt++)
                    standardwave[vdet[ihit].det][vdet[ihit].id]->SetPoint( ipnt, ipnt, standardwave[vdet[ihit].det][vdet[ihit].id]->GetPointY(ipnt) + waveshrink->Eval(ipnt) );
                
                standardwave[vdet[ihit].det][vdet[ihit].id]->SetUniqueID( standardwave[vdet[ihit].det][vdet[ihit].id]->GetUniqueID() + 1 );
            }

        }
        
        if ( vdet.size() > 0 || vref.size() > 0 )
            opt->Fill();

        // display progress and time needed
        if (jentry%100 == 99){
            stop = clock();
            printf("Process %.3f %%  Time remaining %02d min %02d s                                     \r",double(jentry-startentry)/double(stopentry-startentry)*100.,
                int((stop-start)*(stopentry-jentry)/(jentry-startentry)/1e6/60),
                int((stop-start)*(stopentry-jentry)/(jentry-startentry)/1e6)%60);
            fflush(stdout);
        }
    }  // loop for entry
    stop = clock();
    printf("Process %.3f %%  Total Time %02d min %02d s        \n",100.,int((stop-start)/1e6/60),int((stop-start)/1e6)%60);
    
    // TH1/TH2 write
    
    for (int idet = 0; idet < ndet; idet++)
        for (int id = 0; id < nid; id++){
            if ( !valid(idet,id) ) continue;
            
            for (int ipnt = 0; ipnt < 2*length; ipnt++)
                standardwave[idet][id]->SetPoint( ipnt, ipnt, standardwave[idet][id]->GetPointY(ipnt) / standardwave[idet][id]->GetUniqueID() );

            standardwave[idet][id]->SetTitle(Form("det%d id%d nwave=%d", idet, id, standardwave[idet][id]->GetUniqueID()));
            standardwave[idet][id]->Write();
            
            fpol2[idet][id]->Write();
        }
}

void DataAnalysis::BranchOpt()
{
    br_event = 0;

    opt->Branch("ref", &vref);
    opt->Branch("det", &vdet);
}

// 
// event.cc ends here
