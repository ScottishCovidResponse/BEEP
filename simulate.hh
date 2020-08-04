#ifndef BEEPMBP__SIMULATE_HH
#define BEEPMBP__SIMULATE_HH

#include "data.hh"

class DATA;
class MODEL;
class POPTREE;

void simulatedata(DATA &data, MODEL &model, POPTREE &poptree, Mcmc &mcmc);

#endif
