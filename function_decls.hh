void init();                               // Function declarations
void simulatedata();
void PMCMC();
void readdata();
double bootstrap();
void addcomp(string name, double infectivity);
void addparam(string name, double val, double min, double max);
void addtrans(string from, string to, short type, string param1, string param2);
void betaspline();
double sample();
