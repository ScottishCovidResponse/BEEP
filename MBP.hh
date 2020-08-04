#ifndef BEEPMBP__MBP_HH
#define BEEPMBP__MBP_HH

using namespace std;

enum class proposalsmethod
{
	allchainsallparams,
	fixednum,
	fixedtime
};

void MBP(DATA &data, MODEL &model, POPTREE &poptree, Mcmc &mcmc, Mpi &mpi, enum proposalsmethod propmethod);
#endif
