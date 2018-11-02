#include<iostream>
#include<random>
#include<math.h>
#include<vector>
#include<fstream>
#include<tuple>


using namespace std;

double E_i(vector<int>x, vector<int>y){
  double s=0;
  for(int i=0; i<x.size(); i++){  
    s-=2*x[i]*y[i];
    if((i+1)>=x.size()){
      s-=x[i]*x[0];
    }
    else{
      s-=x[i]*x[i+1];
    }
    if((i+1)>=y.size()){
      s-=y[i]*y[0];
    }
    else{
      s-=y[i]*y[i+1];
    }
  }
  return s;
}

double mean_E(double Z, double sum_E){
  return (1/Z)*sum_E;
}

int main(){
  ofstream myfile;
  myfile.open("plot.txt");
  
  int L = 2;
  vector<int> x(L, 0);
  vector<int> y(L, 0);
  vector<int>spin={-1, 1, -1, 1, -1, -1};
  for(int i=0; i<x.size(); i++){
    x[i]=spin[rand()%spin.size()];
    y[i]=spin[rand()%spin.size()];
    //cout << x[i] << "   " <<y[i] << endl;
  }
  double E_b = E_i(x, y);
  double sum_E = E_b*exp(-E_b);
  double Z = exp(-E_b);
  double r;
  cout << mean_E(Z, sum_E)<< endl;
  myfile << 0 << " " << mean_E(Z, sum_E) << " \n";
  
  

for(int i=1; i<1000; i++){
    vector<int> x1 = x;    vector<int> y1 = y;
    int a = rand()%2;
    if(spin[rand()%spin.size()]==1){
      x1[a] = -x1[a];
    }    else{
      y1[a] = -y1[a];
    }
    double E_new = E_i(x1, y1);
    double delta_E = E_new - E_b;
    if(delta_E<=0){
      x=x1;
      y=y1;
      sum_E += E_new*exp(-E_new);
      Z += exp(-E_new);
      E_b = E_new;
      cout << mean_E(Z, sum_E)<< endl;
    }
    else{
      r = (rand()%11)/10;
      double w = exp(-delta_E);
      if(w<=r){
	x=x1;
	y=y1;
	sum_E += E_new*exp(-E_new);
	Z += exp(-E_new);
	E_b = E_new;
	cout << mean_E(Z, sum_E)<< "yyyyyygh" << endl;
      }
    }
    myfile << i << " " << mean_E(Z, sum_E) << " \n";
  }
 double Z_a =12+2*exp(8)+2*exp(-8);
 cout<<"Z:"<<Z_a<<endl;
 double E = (1/Z_a)*(-16*exp(8)+16*exp(-8));
 double E2 = (1/Z_a)*(128*exp(8)+128*exp(-8));
 cout<<"E:"<< E << endl;
 cout<<"Cv:"<< E2-E*E << endl;
 double M2 = (32/Z_a)*(exp(8)+1);
 cout<<"X:"<< M2 << endl;
  return 0;
}
