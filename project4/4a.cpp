#include<iostream>
#include<random>
#include<math.h>
#include<vector>
#include<fstream>


using namespace std;

const int L = 2;

double E_i(int m[L][L]){
  double  s = 0;
  for(int j = 0; j < L; j++){
    for(int i = 0; i < L; i++){
      if((i+2) <= L){
	s += m[j][i]*m[j][i+1];
      }
      else{
	s += m[j][i]*m[j][0];
      }
    }
  }
  for(int i = 0; i < L; i++){
    for(int j = 0; j < L; j++){
      if((j+2) <= L){
	s = m[j][i]*m[j+1][i];
      }
      else{
	s += m[j][i]*m[0][i];
      }
    }
  }
  return s;
}

double mean_E(double Z, double sum_E){
  
  return (1/Z)*sum_E;
}

int main(){ 
  ofstream myfile;
  int m[L][L];
  vector<int>spin={-1, 1, -1, 1, -1, -1};
  for(int i=0; i<L; i++){
    for(int j=0; j<L; j++){
	m[i][j]=spin[rand()%spin.size()];
        
     }

  }
  double E_b = E_i(m);
  double sum_E = E_b*exp(-E_b);
  double Z = exp(-E_b);
  double r;
  cout << mean_E(Z, sum_E)<< endl;
  myfile << 0 << " " << mean_E(Z, sum_E) << " \n";
  
  
  
for(int i=1; i<1000; i++){
    int (*m1)[L] = m;
    int a = rand()%L;
    int b = rand()%L;
   
    m1[b][a] = -m1[b][a];
   
    
    double E_new = E_i(m1);
    double delta_E = E_new - E_b;
    if(delta_E<=0){
       int (*m)[L]  = m1;
       sum_E += E_new*exp(-E_new);
       Z += exp(-E_new);
       E_b = E_new;
       cout << mean_E(Z, sum_E)<< endl;
    }      
      else{
	r = (rand()%11)/10;
	double w = exp(-delta_E);
	if(w<=r){
	     int (*m)[L]  = m1 ;
             sum_E += E_new*exp(-E_new);
	     Z += exp(-E_new);
	     E_b = E_new;
	     cout << mean_E(Z, sum_E)<< "yyyyyygh" << endl;
      }
      } 
      
      myfile << i << " " << mean_E(Z, sum_E) << " \n";	
      
      
  }  
  return 0;
}
