#pragma once

#include <string>

using namespace std;

const int MODE_SIM=0, MODE_PMCMC=1, MODE_MBP=2;         // Different modes of operation 
 
const int MOD_IRISH = 0, MOD_OLD = 1;                   // Different compartmental models

const int FEV_EV=0, INF_EV=1, SET_EV=2;                 // Characterise event types (future event/infection/settime/external)
const int EXT_EV=3, XIFEV_EV=4;

const int EXP_DIST=0, GAMMA_DIST=1;                     // Denotes exponential or gamma distributed 
const int LOGNORM_DIST = 2, INFECTION = 3;

const double tiny = 0.00000001;                         // Used to represent a tiny number
const double large = 1000000;                           // Used to represent a big number

const int checkon = 0;                                  // Set to one to check algorithm is performing correctly
const double finegridsize = 1;//0.02;          // ZZ The range in distance over which the fine grid is used 

const double scale = 684;                  // The number of kilometers across Scotland
const double a = 4.0/scale;                //  Parameters used for spatial kernal
//const double a = 8.0/scale;                //  Parameters used for spatial kernal
const double b = 3;  
const double ddmax = 30.0/scale;         // ZZ The maximum range for kernal
const double rden = 1;                     // Finds the density of houses

const int partmax = 10000;                // The maximum number of particles per core (arbitrarily set)
const int chainmax = 10000;                // The maximum number of chains per core (arbitrarily set)
const int nsettime = 80;                 // The number of time divisions used to represent changes in beta
const int nfix = 1;                      // The number of fixed effects   

const int MAX_NUMBERS = 20000000;

const int SENDRECMAX = 10;        
const int BUFMAX = 2000000;
const int drawfirst = 10;
const double varfac = 4;

const int INFMAX = 50000;                   // The maximum number of infections allowed


