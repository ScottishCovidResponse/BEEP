// Stores details of the simulation / inference

using namespace std;

#include "details.hh"


Mpi::Mpi()
{
	#ifdef USE_MPI
	int num;
	MPI_Comm_size(MPI_COMM_WORLD,&num); ncore = (unsigned int) num;
  MPI_Comm_rank(MPI_COMM_WORLD,&num); core = (unsigned int) num;
	#endif
	
	#ifndef USE_MPI
	ncore = 1;
	core = 0;
	#endif
}

Details::Details(Inputs &inputs)
{
	mode = inputs.mode();          // Gets the mode of operation
	
}
