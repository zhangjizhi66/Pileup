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

#define CHECKFILEPATH "../../CSNSwave202310/"   //Path to the check ROOT file
#define CHECKFILENAME "check"                   //The file name of the check file

#define ROOTFILEPATH "../../CSNSwave202310/"  //The path to generate the ROOT file
#define ROOTFILENAME "standardwave"

const int ndet = 16;
const int nid = 16;

const int L = 2;         // fastfilter parameter
const int nbase = 1000;     // baseline length at the head of waveform
const int trise = 10;    // wave rise time
    
const int length = 100;   // standard waveform length = 2*length
const int pileupzone = 8000;
    
const int minheight = 3000;    // minimum height to build standard wave
const int maxheight = 10000;  // maximum height to build standard wave

const int maxoverflow = 20000;
const int refbin = 10000;

const double standardfwhm = 20;

#define WAVEFORM
// #define ENERGYSUM
// #define QDCSUM
// #define EXTERNALTS


#define EVENTTIMEWINDOWSWIDTH  1000   //ns



#endif /* _USERDEFINE_H_ */

// 
// UserDefine.h ends here
