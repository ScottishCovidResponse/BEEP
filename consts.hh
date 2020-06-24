#pragma once

#include <string>

using namespace std;

const unsigned int MODE_SIM=0, MODE_PMCMC=1, MODE_MBP=2;         // Different modes of operation 
 
const unsigned int MOD_IRISH = 0;                                // Different compartmental models

const unsigned int FEV_EV=0, INF_EV=1;                           // Event types (future event/infection/settime/external)
const unsigned int SET_EV=2, EXT_EV=3, XIFEV_EV=4, XPFEV_EV=5;

const unsigned int NO_DIST=0, EXP_DIST=1, GAMMA_DIST=2;          // Denotes exponential or gamma distributed 
const unsigned int LOGNORM_DIST=3, INFECTION=4;

const double tiny = 0.00000001;                                  // Used to represent a tiny number
const double small = 0.00001;
const double large = 1000000;                                    // Used to represent a big number
const unsigned int UNSET = 999999999;                            // A large unsigned integer to represent "Unset"

const unsigned int checkon = 0;                                  // Set to one to check algorithm is performing correctly

const unsigned int partmax = 10000;                              // The maximum number of particles per core (arbitrarily set)
const unsigned int chainmax = 10000;                             // The maximum number of chains per core (arbitrarily set)

const unsigned int MAX_NUMBERS = 20000000;                       // The maximum buffer size for Send Recv MPI messages
const unsigned int BUFMAX = 2000000;                             // The maximum buffer size for SendI RecvI MPI messages

const double varfac = 4;                                         // A factor which relaxes the observation model

const unsigned int INFMAX = 200000;                              // The maximum number of infections allowed

const unsigned int BOTH=0, PONLY=1, NOT=2;                       // Use to classify particles in MBPs
