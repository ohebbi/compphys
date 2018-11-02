#include<iostream>
#include<random>
#include<math.h>
#include<vector>
#include<fstream>


using namespace std;

const int L = 5;

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
  mean.push_back((1/(pow(temp,2))*((1/Z)*sum_E_heatcap-pow((1/Z)*sum_E,2)))); // C_v : mean[1]
  
mean.push_back((1/Z)*sum_M); // mean[2]
 mean.push_back((1/temp)*((((1/Z)*sum_E_suscept)-pow(mean[2],2)))); // susceptibility : mean[3]
  return mean;
}


int main(){
  ofstream myfile;
  myfile.open("plot.txt");
  int m[L][L];
  vector<int>spin={1, 1, 1, 1, 1, -1};

  for(int i=0; i<L; i++){
    for(int j=0; j<L; j++){
	     m[i][j]=spin[rand()%spin.size()];
    }
  }

  double temp = 1.0; //temperature

  double E_b = -E_i(m);
  double sum_E = E_b*exp(-E_b/temp);
  double sum_M = M_i(m)*exp(-E_b/temp);
  double sum_E_heatcap = pow(E_b,2)*exp(-E_b/temp);
  double sum_E_suscept = pow(M_i(m),2)*exp(-E_b/temp);

  double Z = exp(-E_b/temp);
  double r;
  vector<double> mean = mean_E_M(Z, sum_E, sum_M, sum_E_heatcap, sum_E_suscept, temp);
    myfile << 0 << " " << mean[0] << " " << mean[1] << " " << mean[2] << " " << mean[3] <<  " \n";


for(int i=1; i<1e6; i++){

    int (*m1)[L] = m;
    int a = rand()%L;
    int b = rand()%L;

    m1[b][a] = -m1[b][a];


    double E_new = -E_i(m1);
    double delta_E = E_new - E_b;
    if(delta_E<=0){

       int (*m)[L] = m1;
       sum_E += E_new*exp(-E_new/temp);
       Z     += exp(-E_new/temp);
       sum_M += M_i(m)*exp(-E_new/temp);
       sum_E_heatcap += pow(E_new,2) *exp(-E_new/temp);
       sum_E_suscept += pow(M_i(m1),2)*exp(-E_new/temp);
       E_b    = E_new;
       mean   = mean_E_M(Z, sum_E, sum_M, sum_E_heatcap, sum_E_suscept, temp);

    }
    else{
	     r = (rand()%11)/10.0;
	     double w = exp(-delta_E);

	     if(w<=r){
	        int (*m)[L]  = m1 ;
		sum_E += E_new*exp(-E_new/temp);
	        Z     += exp(-E_new/temp);
	        sum_M += M_i(m)*exp(-E_new/temp);
		sum_E_heatcap += pow(E_new,2) *exp(-E_new/temp);
		sum_E_suscept += pow(M_i(m1),2)*exp(-E_new/temp);
	        E_b    = E_new;
	        mean   = mean_E_M(Z, sum_E, sum_M, sum_E_heatcap, sum_E_suscept, temp);
       }
    }

    myfile << i << " " << mean[0] << " " << mean[1] << " " << mean[2] << " " << mean[3] << " \n";


 }
 cout << mean[0] << " " << mean[1] << " " << mean[2] << " " << mean[3] << " \n";
 myfile.close();
 return 0;
}
