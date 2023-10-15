// UserDefine.h --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 日 1月 10 21:11:47 2021 (+0800)
// Last-Updated: 四 4月 27 10:44:41 2023 (+0800)
//           By: CIAE DAQ
//     Update #: 13
// URL: http://wuhongyi.cn

#ifndef _USERDEFINE_H_
#define _USERDEFINE_H_

#define RAWFILEPATH "../../CSNSwave202310/"   //Path to the original ROOT file
#define RAWFILENAME "event"                   //The file name of the original file

#define WAVEFILEPATH "../../CSNSwave202310/"   //Path to the standardwave ROOT file
#define WAVEFILENAME "standardwave"            //The file name of the standardwave file

#define ROOTFILEPATH "../../CSNSwave202310/"  //The path to generate the ROOT file
#define ROOTFILENAME "ana"

const int ndet = 16;
const int nid = 16;

const int L = 3;         // fastfilter parameter
const int nbase = 1500;     // baseline length at the head of waveform
const int trise = 10;    // wave rise time

// standard waveform range use
const int rangeuseleft = -25;
const int rangeuseright = 75;

const int length = 100;   // standard waveform length = 2*length
const int minpileup = rangeuseright - rangeuseleft;  // pile up peaks will be fitted together

const double overflow = 10000;  // find overflow value at (overflow, +inf)

const int peak2trig = trise;

const int flushlowedge = 2400;

const double standardfwhm = 20;

// fitting parameter
const int maxpeakfitinterval = 2;      // max peak initial value interval
const int targetchi2ndf = 100;         // target chi2/ndf of fitting
const int maxaddpeaks = 5;             // max npeaks allowed being added to fitting
const int lowestfittedpeak = 50;       // if lower, abandon this peak
const double optimizechi2ndf_1 = 0.9;  // add a peak must lower chi2ndf by this portion, otherwise no extra peak
const double optimizechi2ndf_2 = 0.7;  // add a peak must lower chi2ndf by this portion, otherwise no extra peak
const int maxoverflow = 20000;       // upper limit of peak height
const int lowestheight = 5;           // lowest height, times of noise level

#define WAVEFORM
// #define ENERGYSUM
// #define QDCSUM
// #define EXTERNALTS


#define EVENTTIMEWINDOWSWIDTH  1000   //ns



#endif /* _USERDEFINE_H_ */

// 
// UserDefine.h ends here
