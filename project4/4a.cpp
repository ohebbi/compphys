#include<iostream>
#include<random>
#include<math.h>
#include<vector>
#include<fstream>
    

using namespace std;
random_device rd;
mt19937_64 gen(rd());
const int L = 2;
const int nspins = L*L;
uniform_real_distribution<double> dist(0.0, nspins);

const double temp = 1.0; //temperature



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

int spin(double ixx){
  if(ixx<0.5*nspins){
    return 1;
  }
  else{
    return -1;
  }
  }

int main(){
  int high_tech=0;
  int ix;
  ofstream myfile;
  myfile.open("plot.txt");
  int m[L][L];
  for(int i=0; i<L; i++){
    for(int j=0; j<L; j++){
      ix = (double) (dist(gen));
      //cout << spin(ix) << endl;
      m[i][j]=spin(ix);
      
    }
  }

  
  double E_b = -E_i(m);
  double B = exp(-E_b/temp);
  double sum_E = E_b*B;
  double sum_M = M_i(m)*B;
  double sum_E_heatcap = pow(E_b,2)*B;
  double sum_E_suscept = pow(M_i(m),2)*B;
  int teller=0;

  double Z = exp(-E_b/temp);
  
  double r;
  vector<double> mean = mean_E_M(Z, sum_E, sum_M, sum_E_heatcap, sum_E_suscept, temp);
  myfile << 0 << " " << mean[0] << " " << mean[1] << " " << mean[2] << " " << mean[3] <<  " \n"; 


  for(int i=1; i<1e6; i++){

    int (*m1)[L] = m;
    int a = (int) (dist(gen)/L);
    int b = (int) (dist(gen)/L);

    m1[b][a] = -m1[b][a];

    double E_new = -E_i(m1);
    double delta_E = E_new - E_b;
    if(delta_E<=0){
      teller+=1;
       B=exp(-E_new/temp);
       int (*m)[L] = m1;
       sum_E += E_new*B;
       Z     += B;
       sum_M += M_i(m)*B;
       sum_E_heatcap += pow(E_new,2) *B;
       sum_E_suscept += pow(M_i(m1),2)*B;
       E_b    = E_new;
       mean   = mean_E_M(Z, sum_E, sum_M, sum_E_heatcap, sum_E_suscept, temp);

    }
    else{
       r = ( (double) (dist(gen))/nspins);
       cout << r << endl;
       double w = exp(-delta_E);
       B=exp(-E_new/temp);
       if(w<=r){
	 teller+=1;
	 int (*m)[L]  = m1 ;
	 sum_E += E_new*B;
	 Z     += B;
	 sum_M += M_i(m)*B;
	 sum_E_heatcap += pow(E_new,2) *B;
       	 sum_E_suscept += pow(M_i(m1),2)*B;
	 E_b    = E_new;
	 mean   = mean_E_M(Z, sum_E, sum_M, sum_E_heatcap, sum_E_suscept, temp);
       }
    }
    if(i%1==0){
      myfile << i << " " << mean[0] << " " << mean[1] <<" "<< mean[2] << " "<< mean[3] << " \n";
    }
    if(i>6e8){
      myfile << mean[1] <<" \n";
    }
    if(i>1e8 and i%10000000==0){
      cout << high_tech << "%"<< endl;
      high_tech+=1;
    }

 }
  cout << mean[0] << " " << mean[1] << " "<< teller << " \n";
 myfile.close();
 return 0;
}
