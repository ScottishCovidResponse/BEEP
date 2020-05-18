#pragma once

void definemodel();
void betaspline(long nsettime, short tmax, const vector <double> &splinet, short nspline, const vector <PARAM> &param,
								double settime[nsettime], double beta[nsettime]);
