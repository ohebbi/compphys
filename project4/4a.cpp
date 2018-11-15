#include<iostream>
#include<random>
#include<math.h>
#include<vector>
#include<fstream>
    

using namespace std;


const int L = 40;
const int nspins = L*L;
const int n = 1.0e8;
const double temp = 1.0; //temperature


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

int main(){
ofstream myfile;
myfile.open("plot.txt");
//for(double temp = 2.0; temp <= 2.3; temp += 0.05){
  double Eaverage;
  double E2average;
  double Maverage;
  double M2average;
  double Mabsaverage;
  double ix;

  vector<double> prob;
 
  double w[17];
  double mean[5]; 
  

  
  int m[L][L];
  for(int i=0; i<L; i++){
    for(int j=0; j<L; j++){
      ix = dis(gen);          
      m[i][j]=spin(ix);
      cout << m[i][j] << endl;
     
    }
  }

  for(int dE = -8; dE <= 8; dE++) w[dE+8] = 0;
  for(int dE = -8; dE <= 8; dE+=4) w[dE+8] = exp(-dE/temp);
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
myfile << 1 << " " << Eaverage/nspins  << " " << (E2average-Eaverage*Eaverage)/(temp*temp*nspins) << " " << Mabsaverage/nspins << " " << (M2average-Maverage*Maverage)/(temp*nspins) <<   " \n";
     
   
  for(int i=2; i<n; i++){

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
  myfile << i << " " << Eaverage/nspins  << " " << (E2average-Eaverage*Eaverage)/(nspins) << " " << Mabsaverage/nspins << " " << (M2average-Maverage*Maverage)/(temp*nspins) <<  "\n";
}
}
/*
Eaverage = mean[0]/n;
E2average = mean[1]/n;
Maverage = mean[2]/n;
M2average = mean[3]/n;
Mabsaverage = mean[4]/n;
myfile << temp << " " << Eaverage/nspins  << " " << (E2average-Eaverage*Eaverage)/(temp*temp*nspins) << " " << Mabsaverage/nspins << " " << (M2average-Maverage*Mabsaverage)/(temp*nspins) << " \n";
}*/

 myfile.close();
 return 0;
}
