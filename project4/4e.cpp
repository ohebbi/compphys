#include<iostream>
#include<random>
#include<math.h>
#include<vector>
#include<fstream>
#include <mpi.h>


using namespace std;

double E;
double T_init=2.0;
double T_final=2.3;
double delta_T=0.01;
double L=20;

double mean_E(double Ti){
  return E=Ti; 
}

int main(int nargs, char* args[]){
  int numprocs, my_rank; 
  double  time_start, time_end, total_time;
  
  ofstream myfile;
  myfile.open("plot1.txt");
  
  MPI_Init (&nargs, &args);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  time_start = MPI_Wtime();
  // all the calculations in here
  double interval = (T_final-T_init)/numprocs;
  double T_0=T_init + my_rank*interval;
  double T_1=T_final-(numprocs-my_rank-1)*interval;
  double n = (T_final-T_init)/(delta_T*numprocs);

  for(int i=0; i<=n; i++){
    myfile << T_0 << "  "<< mean_E(T_0)  <<  " \n";
    T_0+=delta_T;
  }
 
  time_end = MPI_Wtime();
  total_time = time_end-time_start;

  MPI_Finalize (); 
  return 0;
  
}

