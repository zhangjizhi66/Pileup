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
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TSpectrum.h"
#include "TMarker.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <map>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TGraph *standardwave[ndet][nid];
TF1 *fpol2[ndet][nid];

int fitdet, fitid;
double rejval;
std::vector<int> vrejpnts;

double ffit(double *val, double *par)
{
    double x0 = val[0];
    int npeaks = par[0];
    
    for ( int rejpnt : vrejpnts )
        if ( abs(x0 - rejpnt) < 1 ){
            TF1::RejectPoint();
            return rejval;
       }
    
    double wave = 0;
    for (int ipeak = 0; ipeak < npeaks; ipeak++){
        double A = par[3*ipeak+1];
        double pos = par[3*ipeak+2];
        //int type = par[3*ipeak+3];
        
        double wavefwhm = fpol2[fitdet][fitid]->Eval(A);
        
        double x = pos + (x0-pos)*standardfwhm/wavefwhm;
        if ( x - pos < rangeuseleft || x - pos > rangeuseright )
            wave += 0;
        else
            wave += A * standardwave[fitdet][fitid]->Eval(x - pos + length);
    }
    return wave;
}

void DataAnalysis::Init()
{
    if (ipt == NULL) return;
    ipt->SetBranchAddress("event", &br_event);
  
    // init wave
    for (int idet = 0; idet < ndet; idet++)
        for (int id = 0; id < nid; id++){
            if ( !valid(idet,id) ) continue;
            
            standardwave[idet][id] = (TGraph*)ipf_wave->Get(Form("standardwave_det%d_id%d",idet,id));
            fpol2[idet][id] = (TF1*)ipf_wave->Get(Form("fpol2_det%d_id%d",idet,id));
        }
}

void DataAnalysis::Loop(TFile *opf_, TTree *opt_, Long64_t startentry = -1, Long64_t stopentry = -1)
{
    if (opt_ == NULL) return;
    
    opf_->cd();

    opt = opt_;
    BranchOpt();
    
    clock_t start = clock(), stop = clock();
    
    if (startentry == -1 && stopentry == -1){
        startentry = 0;
        stopentry = ipt->GetEntries();
    }
    
    if ( stopentry > ipt->GetEntries() )
        stopentry = ipt->GetEntries();
  
    for (Long64_t jentry = startentry; jentry < stopentry; jentry++){
        ipt->GetEntry(jentry);
        // if (jentry > 1000) continue;
       
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
            det.twave.clear();
            det.ewave.clear();
            det.nwave.clear();
#endif
            //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

            if ( flag == 0 ){     // ref for beam bursts
                vref.clear();
                vref.push_back(det);
                //continue;
            }
            else
                if ( (*br_event)[i].ltra > 0 ) vdet.push_back(det);
        }
        
        for (int ihit = 0; ihit < int(vdet.size()); ihit++){
            if ( int(vdet[ihit].wave.size()) <= 0 ) continue;
            
            fitdet = vdet[ihit].det;
            fitid = vdet[ihit].id;
            
            if ( !valid(fitdet,fitid) ) continue;
            
            // get baseline and noise
            
            double base = 0;
            for (int ipnt = 0; ipnt < nbase; ipnt++)
                base += (double)vdet[ihit].wave[ipnt] / nbase;
            
            double noise = 0;
            for (int ipnt = 0; ipnt < nbase; ipnt++)
                if ( abs( vdet[ihit].wave[ipnt] - base ) > noise )
                    noise = abs( vdet[ihit].wave[ipnt] - base );
            
            // store as TGraph
            TGraph *gwave = new TGraph;
            for (int ipnt = 0; ipnt < int(vdet[ihit].wave.size()); ipnt++)
                gwave->SetPoint(ipnt, ipnt, vdet[ihit].wave[ipnt] - base);
            
            // find overflow value and reject points
            rejval = -10000;
            for (int ipnt = 1; ipnt < gwave->GetN()-1; ipnt++)
                if ( gwave->GetPointY(ipnt) == gwave->GetPointY(ipnt-1) && gwave->GetPointY(ipnt) == gwave->GetPointY(ipnt+1) && gwave->GetPointY(ipnt) > overflow ){
                    rejval = gwave->GetPointY(ipnt);
                    break;
                }
            vrejpnts.clear();
            for ( int ipnt = 0; ipnt < gwave->GetN(); ipnt++ )
                if ( gwave->GetPointY(ipnt) == rejval )
                    vrejpnts.push_back(ipnt);
            
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
            
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            
            
            int npeaks0 = 1;
            for (int itrig = 0; itrig < int(vtrig.size()); itrig += npeaks0){

                if ( vtrig[itrig] < flushlowedge ) continue;
                
                // get number of peaks (pileup)
                npeaks0 = 1;
                for (int jtrig = itrig+1; jtrig < int(vtrig.size()); jtrig++){
                    if (vtrig[jtrig] - vtrig[jtrig-1] <= minpileup)
                        npeaks0++;
                    else
                        break;
                }

                // fit range
                int left = std::max(0, vtrig[itrig+0] + rangeuseleft);
                int right = std::min( vtrig[itrig+npeaks0-1] + rangeuseright, int(vdet[ihit].wave.size()) );
                right += 10;
                
                // lowest height cut
                double height = 0;
                for ( int ipnt = left; ipnt < right; ipnt++ )
                    if ( gwave->GetPointY(ipnt) > height )
                        height = gwave->GetPointY(ipnt);
                
                if ( height < lowestheight*noise ) continue;
                
                // find peaks
                std::vector<double> vx, vy;
                for (int jtrig = itrig; jtrig < itrig+npeaks0; jtrig++)
                    vx.push_back( vtrig[jtrig] + peak2trig );
                
                // fit
                TF1 *f = new TF1("ffit", ffit, left, right, 3*npeaks0+1);
                f->SetNpx(right-left);
                double chi2ndf = -1, chi2ndf1 = -1;
                TFitResultPtr fr;
                std::multimap<Int_t,Double_t> me;
                int npeaks = vx.size();
                
                // first fit
                f->FixParameter(0, npeaks);
                for (int ipeak = 0; ipeak < npeaks; ipeak++){
                    f->SetParameter(3*ipeak+1, gwave->GetPointY(vx[ipeak]));
                    f->SetParLimits(3*ipeak+1, 0, 100*rejval);
                    f->SetParameter(3*ipeak+2, vx[ipeak]);
                    f->SetParLimits(3*ipeak+2, left - rangeuseright, right - rangeuseleft);
                    f->FixParameter(3*ipeak+3, 0);
                }
                
                fr = gwave->Fit(f, "SQRN", "", left, right);
                chi2ndf = fr->Chi2() / fr->Ndf();
                
                for (int ipeak = 0; ipeak < npeaks; ipeak++){
                    vy.push_back(f->GetParameter(3*ipeak+1));
                    vx[ipeak] = f->GetParameter(3*ipeak+2);
                }
                
                do {
                    npeaks = vx.size();
                    
                    // find residual error peaks
                    TH1I* hdiff = new TH1I("", "", right-left, left, right);
                    for (int ipnt = left; ipnt < right; ipnt++)
                        hdiff->SetBinContent(ipnt-left, abs( gwave->GetPointY(ipnt) - f->Eval(ipnt) ) );
                    TSpectrum *s = new TSpectrum(500);
                    int nfound = s->Search(hdiff, 2, "", 0.1);
                    double *xpeaks = s->GetPositionX();
                    double *ypeaks = s->GetPositionY();
                    
                    me.clear();
                    for (int j = 0; j < nfound; j++){
                        if ( ypeaks[j] < lowestheight*noise ) continue;
                        
                        bool repeat = false;
                        for (int ipeak = 0; ipeak < npeaks; ipeak++)
                            if ( abs( xpeaks[j] - vx[ipeak] ) <= maxpeakfitinterval ){
                                repeat = true;
                                break;
                            }
                        if (repeat) continue;
                        
                        me.insert(std::make_pair(int(ypeaks[j]),xpeaks[j]));
                    }
                    if ( me.empty() ) break;
                    
                    // add new peak and fit
                    bool addnewpeak = false;
                    if (f) delete f;
                    f = new TF1("ffit", ffit, left, right, 3*(npeaks+1)+1);
                    f->SetNpx(right-left);
                    f->FixParameter(0, npeaks+1);
                    for (auto im = me.rbegin(); im != me.rend(); im++){
                        for (int ipeak = 0; ipeak < npeaks; ipeak++){
                            if ( vx[ipeak] - im->second <= rangeuseleft || vx[ipeak] - im->second >= rangeuseright ){
                                f->FixParameter(3*ipeak+1, vy[ipeak]);
                                f->FixParameter(3*ipeak+2, vx[ipeak]);
                            }
                            else {
                                f->SetParameter(3*ipeak+1, vy[ipeak]);
                                f->SetParLimits(3*ipeak+1, 0, 100*rejval);
                                f->SetParameter(3*ipeak+2, vx[ipeak]);
                                f->SetParLimits(3*ipeak+2, left - rangeuseright, right - rangeuseleft);
                            }
                            f->FixParameter(3*ipeak+3, 0);
                        }
                        f->SetParameter(3*npeaks+1, gwave->GetPointY(im->second));
                        f->SetParLimits(3*npeaks+1, 0, 100*rejval);
                        f->SetParameter(3*npeaks+2, im->second);
                        f->SetParLimits(3*npeaks+2, left - rangeuseright, right - rangeuseleft);
                        f->FixParameter(3*npeaks+3, 0);

                        fr = gwave->Fit(f, "SQRN", "", left, right);
                        chi2ndf1 = fr->Chi2() / fr->Ndf();
                        
                        if ( chi2ndf1 < optimizechi2ndf*chi2ndf ){
                            addnewpeak = true;
                            break;
                        }
                    }
                    
                    if ( addnewpeak ){
                        for (int ipeak = 0; ipeak < npeaks; ipeak++){
                            vx[ipeak] = f->GetParameter(3*ipeak+2);
                            vy[ipeak] = f->GetParameter(3*ipeak+1);
                        }
                        vx.push_back(f->GetParameter(3*npeaks+2));
                        vy.push_back(f->GetParameter(3*npeaks+1));
                        chi2ndf = chi2ndf1;
                    }
                    else
                        break;
                } while (chi2ndf > targetchi2ndf && int(vx.size()) <= npeaks0 + maxaddpeaks);
                
                bool clean = false;
                while ( !clean ){
                    int ip = 0;
                    if (f) delete f;
                    f = new TF1("ffit", ffit, left, right, 3*vx.size()+1);
                    f->SetNpx(right-left);
                    clean = true;
                    for (int ipeak = 0; ipeak < int(vx.size()); ipeak++){
                        bool repeat = false;
                        for (int jp = 0; jp < f->GetParameter(0); jp++)
                            if ( abs(vx[ipeak] - f->GetParameter(3*jp+2)) <= maxpeakfitinterval )
                                repeat = true;
                        if ( vy[ipeak] > lowestfittedpeak && !repeat ){
                            f->SetParameter(3*ip+1, vy[ipeak]);
                            f->SetParLimits(3*ip+1, 0, 100*rejval);
                            f->SetParameter(3*ip+2, vx[ipeak]);
                            f->SetParLimits(3*ip+2, vx[ipeak]-maxpeakfitinterval, vx[ipeak]+maxpeakfitinterval);
                            f->FixParameter(3*ip+3, 0);
                            f->FixParameter(0, ++ip);
                        }
                        else
                            clean = false;
                    }
                    for (int ipar = 3*ip+1; ipar < 3*int(vx.size())+1; ipar++)
                        f->FixParameter(ipar, 0);

                    if (clean)
                        fr = gwave->Fit(f, "SQR+", "", left, right);
                    else
                        fr = gwave->Fit(f, "SQRN", "", left, right);
                    chi2ndf = fr->Chi2() / fr->Ndf();
                    
                    vx.clear();
                    vy.clear();
                    for (int ipeak = 0; ipeak < f->GetParameter(0); ipeak++){
                        vy.push_back( f->GetParameter(3*ipeak+1) );
                        vx.push_back( f->GetParameter(3*ipeak+2) );
                    }
                }
                
                for (int ipeak = 0; ipeak < f->GetParameter(0); ipeak++){
                    vdet[ihit].ewave.push_back( f->GetParameter(3*ipeak+1) );
                    vdet[ihit].twave.push_back( f->GetParameter(3*ipeak+2) );
                    vdet[ihit].nwave.push_back( f->GetParameter(0) );
                    
                    TMarker *mark = new TMarker(f->GetParameter(3*ipeak+2), f->GetParameter(3*ipeak+1), 29);
                    mark->SetMarkerColor(kRed);
                    mark->SetMarkerSize(1);
                    gwave->GetListOfFunctions()->Add(mark);
                }
                
                if (f) delete f;
            }
            
            if ( jentry < 10 ){
                gwave->SetName(Form("gwave_%d_%d_%lld_%d", vdet[ihit].det, vdet[ihit].id, jentry, ihit));
                gwave->Write();
            }
            delete gwave;
        }
        
        if ( vdet.size() > 0 )
            opt->Fill();

        // display progress and time needed
        if (jentry%2 == 1){
            stop = clock();
            printf("Process %.3f %%  Time remaining %02d min %02d s                                     \r",double(jentry-startentry)/double(stopentry-startentry)*100.,
                int((stop-start)*(stopentry-jentry)/(jentry-startentry)/1e6/60),
                int((stop-start)*(stopentry-jentry)/(jentry-startentry)/1e6)%60);
            fflush(stdout);
        }

    }  // loop for entry

    stop = clock();
    printf("Process %.3f %%  Total Time %02d min %02d s        \n",100.,int((stop-start)/1e6/60),int((stop-start)/1e6)%60);
    
    opt->Write();
    
    // TH1/TH2 write
}

void DataAnalysis::BranchOpt()
{
    br_event = 0;

    opt->Branch("ref", &vref);
    opt->Branch("det", &vdet);
}

// 
// event.cc ends here
