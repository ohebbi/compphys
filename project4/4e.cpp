#include<iostream>
#include<random>
#include<math.h>
#include<vector>
#include<fstream>
#include <mpi.h>


using namespace std;


const int L = 40;
const int nspins = L*L;

const int N = 1e6;


const double T_init=2.0;
const double T_final=2.3;
const double delta_T=0.02;

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<double> dis(0.0, 1.0);
uniform_real_distribution<double> dist(0.0, L);

inline int per_bound(int i, int limit, int add){
  return(i+limit+add) % (limit);
}

double M_i(int m[L][L]){
  double M = 0;
  for(int i = 0; i <L; i++){
    for(int j = 0; j < L; j++){
      M += m[i][j];
    }
  }
  return M;
}

double E_i(int m[L][L]){
 double  s = 0;

  for(int j = 0; j < L; j++){
    for(int i = 0; i < L; i++){
      if((i+1) < L){
	      s += m[j][i]*m[j][i+1];
      }
      else{
	      s += m[j][L-1]*m[j][0];
      }
    }
  }

  for(int i = 0; i < L; i++){
    for(int j = 0; j < L; j++){
      if((j+1) < L){
	      s += m[j][i]*m[j+1][i];
      }
      else{
	      s += m[L-1][i]*m[0][i];
      }
    }
  }
  return s;
}

vector<double> Metropolis(int l, int m[L][L], double E, double M, double *w){
  for(int x = 0; x < l; x++){
    for(int y = 0; y < l; y++){
      int rx = dist(gen);
      int ry = dist(gen);
      int dE = 2*m[rx][ry]*(m[rx][per_bound(ry, l, -1)]+m[per_bound(rx,l,-1)][ry]+m[rx][per_bound(ry,l,1)]+m[per_bound(rx,l,1)][ry]);

      if(dis(gen) <= w[dE+8]){
        m[rx][ry] *= -1;
        M += (double)2*m[rx][ry];
        E += (double)dE;

      }

    }
  }
  return {E, M};
}



int spin(double ixx){
  if(ixx>=0.5){
    return 1;
  }
  else{
    return -1;
  }
}

int main(int nargs, char* args[]){

  double ix;
  int a;
  int b;
  int t;
  double E_new;
  double Eaverage;
  double E2average;
  double Maverage;
  double M2average;
  double Mabsaverage;

  vector<double> prob;

  double w[17];
  double mean[5];


  int numprocs, my_rank;
  double  time_start, time_end, total_time;

  ofstream myfile1;
  myfile1.open("L100-1.txt");
  ofstream myfile2;
  myfile2.open("L100-2.txt");

  int m[L][L];
  for(int i=0; i<L; i++){
    for(int j=0; j<L; j++){
      ix = dis(gen);
      m[i][j]=spin(ix);
    }
  }

  MPI_Init (&nargs, &args);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  time_start = MPI_Wtime();

  double interval = (T_final-T_init)/numprocs;
  double T_0=T_init+my_rank*interval;
  double T_1=T_final-(numprocs-my_rank-1)*interval;
  double n=(T_final-T_init)/(delta_T*numprocs);

  for(int dE = -8; dE <= 8; dE++) w[dE+8] = 0;
  for(int dE = -8; dE <= 8; dE+=4) w[dE+8] = exp(-dE/T_0);

  float E = -E_i(m);
  double M = M_i(m);

  vector<double> EM;

  mean[0] = E;
  mean[1] = E*E;
  mean[2] = M;
  mean[3] = M*M;
  mean[4] = fabs(M);

  Eaverage = mean[0];
  E2average = mean[1];
  Maverage = mean[2];
  M2average = mean[3];
  Mabsaverage = mean[4];
  if(my_rank==0){
  myfile1 << 1 << " " << Eaverage/nspins  << " " << (E2average-Eaverage*Eaverage)/(T_0*T_0*nspins) << " " << Mabsaverage/nspins << " " << (M2average-Maverage*Maverage)/(T_0*nspins) <<   " \n";
  }

  for(int j=0; j<=n; j++){

  for(int i=2; i<N; i++){

    EM = Metropolis(L, m, E, M, w);
    float Echeck = E;
    E = EM[0];

    M = EM[1];

    mean[0] += E;
    mean[1] += E*E;
    mean[2] += M;
    mean[3] += M*M;
    mean[4] += fabs(M);


if (i%10000 == 0){
  Eaverage = mean[0]/i;
  E2average = mean[1]/i;
  Maverage = mean[2]/i;
  M2average = mean[3]/i;
  Mabsaverage = mean[4]/i;
  if(my_rank==0){
   myfile1 << T_0  << " " << Eaverage/nspins  << " " << (E2average-Eaverage*Eaverage)/(T_0*T_0*nspins) << " " << Mabsaverage/nspins << " " << (M2average-Maverage*Maverage)/(T_0*nspins) <<   " \n";

  }
  else{
    myfile2 << T_0  << " " << Eaverage/nspins  << " " << (E2average-Eaverage*Eaverage)/(T_0*T_0*nspins) << " " << Mabsaverage/nspins << " " << (M2average-Maverage*Maverage)/(T_0*nspins) <<   " \n";
  }
}
}
  //cout << T_0 << mean[0] << " " << mean[1] << i << " " << mean[2] << " " << mean[3] << my_rank << " \n";

  cout<< T_0<<"\n";
  T_0+=delta_T;

 }
 //cout << T_0 << mean[0] << " " << mean[1] << i << " " << mean[2] << " " << mean[3] << my_rank << " \n";
  myfile1.close();
  myfile2.close();

  time_end = MPI_Wtime();
  total_time = time_end-time_start;

  MPI_Finalize ();
  cout << '\a';
  return 0;

}
