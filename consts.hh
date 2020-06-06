#pragma once

using namespace std;

const string type = "small";               // These are different models that have been analysed
//const string type = "med";
//const string type = "realscotsim";
//const string type = "scotsim";
//const string type = "uksim";
 
const short MOD_IRISH = 0, MOD_OLD = 1;    // Different compartmental models
const short modelsel = MOD_IRISH;

const short FEV_EV=0, INF_EV=1, SET_EV=2;  // Used to characterise an event type (future event/infection/settime/external)
const short EXT_EV=3, XIFEV_EV=4;
const double tiny = 0.00000001;            // Used to represent a tiny number
const double large = 1000000;              // Used to represent a big number
const short timestep = 7;                  // Case data uses week interval

const short checkon = 0;                   // Set to one to check algorithm is performing correctly
const double finegridsize = 1;//0.02;          // ZZ The range in distance over which the fine grid is used 

const double scale = 684;                  // The number of kilometers across Scotland
//const double a = 4.0/scale;                //  Parameters used for spatial kernal
const double a = 8.0/scale;                //  Parameters used for spatial kernal
const double b = 2;     // ZZ
const double ddmax = 30.0/scale;         // ZZ The maximum range for kernal
const double rden = 1;                     // Finds the density of houses

const long fediv=17*7*10;                   // The number of divisions into which the global timeline is divided

const long partmax = 10000;                // The maximum number of particles (arbitrarily set)
const short EXP_DIST=0, GAMMA_DIST=1;      // Denotes exponential or gamma distributed 
const short LOGNORM_DIST = 2, INFECTION = 3;
const long nsettime = 80;                 // The number of time divisions used to represent changes in beta
const short nfix = 1;                      // The number of fixed effects   

const int MAX_NUMBERS = 10000000;

const short SENDRECMAX = 10;        
const long BUFMAX = 1000000;
