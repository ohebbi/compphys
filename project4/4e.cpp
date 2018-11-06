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
double delta_T=0.05;

const int L = 20;

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

vector<double> mean_E_M(double Z, double sum_E, double sum_M, double sum_E_heatcap, double sum_E_suscept, double temp){
  vector<double> mean;
  mean.push_back((1/Z)*sum_E); // mean[0]
  //mean.push_back((1/(pow(temp,2))*((1/Z)*sum_E_heatcap-pow((1/Z)*sum_E,2)))); // C_v : mean[1]
  
mean.push_back((1/Z)
vector<double> mean_E_*sum_M); // mean[2]
//mean.push_back((1/temp)*((((1/Z)*sum_E_suscept)-pow(mean[2],2)))); // susceptibility : mean[3]
  return mean;
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

  double interval = (T_final-T_init)/numprocs;
  double T_0=T_init+my_rank*interval;
  double T_1=T_final-(numprocs-my_rank-1)*interval;
  double n=(T_final-T_init)/(delta_T*numprocs);

  int m[L][L];
  vector<int>spin={1, -1, -1, -1, 1, -1};

  for(int i=0; i<=n; i++){

  for(int i=0; i<L; i++){
    for(int j=0; j<L; j++){
	     m[i][j]=spin[rand()%spin.size()];
    }
  }

  double temp = 1.0; //temperature
  
  double E_b = -E_i(m);
  double B = exp(-E_b/temp);
  double sum_E = E_b*B;
  double sum_M = M_i(m)*B;
  double sum_E_heatcap = pow(E_b,2)*B;
  double sum_E_suscept = pow(M_i(m),2)*B;

  double Z = exp(-E_b/temp);
  double r;
  vector<double> mean = mean_E_M(Z, sum_E, sum_M, sum_E_heatcap, sum_E_suscept, temp);
    myfile << 0 << " " << mean[0] << " " << mean[1] << " " << mean[2] << " " << mean[3] <<  " \n";


  for(int i=1; i<1e5; i++){

    int (*m1)[L] = m;
    int a = rand()%L;
    int b = rand()%L;

    m1[b][a] = -m1[b][a];

    double E_new = -E_i(m1);
    double delta_E = E_new - E_b;
    if(delta_E<=0){
       B=exp(-E_new/temp);
       int (*m)[L] = m1;
       sum_E += E_new*B;
       Z     += B;
       sum_M += M_i(m)*B;
       //sum_E_heatcap += pow(E_new,2) *B;
       //sum_E_suscept += pow(M_i(m1),2)*B;
       E_b    = E_new;
       mean   = mean_E_M(Z, sum_E, sum_M, sum_E_heatcap, sum_E_suscept, temp);

    }
    else{
       r = (rand()%11)/10.0;
       double w = exp(-delta_E);
       B=exp(-E_new/temp);
       if(w<=r){
	 int (*m)[L]  = m1 ;
	 sum_E += E_new*B;
	 Z     += B;
	 sum_M += M_i(m)*B;
	 //sum_E_heatcap += pow(E_new,2) *B;
       	 //sum_E_suscept += pow(M_i(m1),2)*B;
	 E_b    = E_new;
	 mean   = mean_E_M(Z, sum_E, sum_M, sum_E_heatcap, sum_E_suscept, temp);
       }
    }

    myfile << T_0 << " " << mean[0] << " " << mean[1] << " " << mean[2] << " " << mean[3] << i <<" \n";
    T_0+=delta_T;
  }

 }
  //  cout << mean[0] << " " << mean[1] << " " << mean[2] << " " << mean[3] << my_rank << " \n";
 myfile.close();


  
  time_end = MPI_Wtime();
  total_time = time_end-time_start;

  MPI_Finalize (); 
  return 0;
  
}

