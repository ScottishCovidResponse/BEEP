#pragma once

using namespace std;

//const string type = "small";
//const string type = "med";
const string type = "scotsim";
//const string type = "uksim";
 
const short FEV_EV=0, INF_EV=1, SET_EV=2;  // Used to characterise an event type (future event/infection/settime/external)
const short EXT_EV=3;
const double tiny = 0.00000001;            // Used to represent a tiny number
const double large = 10000000;             // Used to represent a big number
const short timestep = 7;                  // Case data uses week interval

const short checkon = 0;                   // Set to one to check algorithm is performing correctly
const double finegridsize = 0.02;          // The range in distance over which the fine grid is used 
const double d0 = 0.05;                    // The minumum distance cut-off for the matrix M

const long fediv = 1050;                   // The number of divisions into which the global timeline is divided
const long partmax = 10000;                // The maximum number of particles (arbitrarily set)
const short EXP_DIST=0, GAMMA_DIST=1;      // Denotes exponential or gamma distributed 
const long nsettime = 100;                 // The number of time divisions used to represent changes in beta
const short nfix = 1;                      // The number of fixed effects   
