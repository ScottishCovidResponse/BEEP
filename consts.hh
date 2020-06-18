#pragma once

#include <string>

using namespace std;

const int MODE_SIM=0, MODE_PMCMC=1, MODE_MBP=2;         // Different modes of operation 
 
const int MOD_IRISH = 0;                                // Different compartmental models

const int FEV_EV=0, INF_EV=1, SET_EV=2;                 // Characterise event types (future event/infection/settime/external)
const int EXT_EV=3, XIFEV_EV=4, XPFEV_EV=5;

const int NO_DIST=-1, EXP_DIST=0, GAMMA_DIST=1;         // Denotes exponential or gamma distributed 
const int LOGNORM_DIST = 2, INFECTION = 3;

const double tiny = 0.00000001;                         // Used to represent a tiny number
const double large = 1000000;                           // Used to represent a big number

const int checkon = 0;                                  // Set to one to check algorithm is performing correctly

const int partmax = 10000;                              // The maximum number of particles per core (arbitrarily set)
const int chainmax = 10000;                             // The maximum number of chains per core (arbitrarily set)
const int nsettime = 100;                               // The number of time divisions used to represent changes in beta

const int MAX_NUMBERS = 20000000;                       // The maximum buffer size for Send Recv MPI messages
const int BUFMAX = 2000000;                             // The maximum buffer size for SendI RecvI MPI messages

const double varfac = 4;                                // A factor which relaxes the observation model

const int INFMAX = 50000;                               // The maximum number of infections allowed

const int BOTH = 0, PONLY = 1, NOT = 2;                 // Use to classify particles in MBPs
