#pragma once

const short FEV_EV=0, INF_EV=1, SET_EV=2;  // Used to characterise an event type (future event/infection/settime)
const double tiny = 0.00000001;            // Used to represent a tiny number
const double large = 10000000;             // Used to represent a big number
const double invT = 1;                     // The inverse temperature (used to relax the observation model)
//const long popsize = 5000000;              // Number of individuals in Scotland (approximately)
//const long nhouse = 1500000;               // Number of houses
const long popsize = 10000;                 // A smaller test case
const long nhouse = 3000;               
const short checkon = 0;                   // Set to one to check algorithm is performing correctly
const double finegridsize = 0.02;          // The range in distance over which the fine grid is used 
const double d0 = 0.05;                    // The minumum distance cut-off for the matrix M
const short tmax = 105;                    // The time over which simulation / inference is performed
const short RX = 1, RY = 1;                // When siumlating this gives a hypothetical grid of regions
const short nregion = RX*RY;               // When the real data is analysed these will be Healthboard level regions
const long fediv = 1000;                   // The number of divisions into which the global timeline is divided
const long partmax = 10000;                // The maximum number of particles (arbitrarily set)
const short EXP_DIST=0, GAMMA_DIST=1;      // Denotes exponential or gamma distributed 
const long nsettime = 100;                 // The number of time divisions used to represent changes in beta
