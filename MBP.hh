#ifndef BEEPMBP__MBP_HH
#define BEEPMBP__MBP_HH

using namespace std;

enum class proposalsmethod
{
	allchainsallparams,
	fixednum,
	fixedtime
};

void MBP(DATA &data, MODEL &model, POPTREE &poptree, unsigned int nsamp, unsigned int core, unsigned int ncore, unsigned int npart, enum proposalsmethod propmethod);
#endif
