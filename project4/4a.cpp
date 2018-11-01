#include<iostream>
#include<random>
#include<math.h>
#include<vector>
#include<fstream>


using namespace std;

const int L = 2;

double E_i(int m[L][L]){
  double  s = 0;
  for(int i =0; i < L; i++){
    for(int j = 0; j < L; j++){
      
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
	cout << m[i][j];
     }

  }
  double E_b = E_i(m);
  double sum_E = E_b*exp(-E_b);
  double Z = exp(-E_b);
  double r;
  cout << mean_E(Z, sum_E)<< endl;
  myfile << 0 << " " << mean_E(Z, sum_E) << " \n";
  
  
  
for(int i=1; i<100; i++){
    int (*m1)[L] = m;
    int a = rand()%2;
    if(spin[rand()%spin.size()]==1){
      m1[0][a] = -m1[0][a];
    }   
    else{
      m1[1][a] = -m1[1][a];
    }
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
