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

  //tuple<vector<int>> matrix = ({1,0}, {0,1});
  auto matrix = make_tuple (x);
  //for(int i=0; i<L; i++){
    //get<i>(matrix) = 1;
    //cout<<get<i>(matrix)<<endl;
  //}
  vector<int>spin={-1, 1, -1, 1, -1, -1};
  /*
  for(int i=0; i<x.size(); i++){
    for(int j=0; i<matrix.size(); j++){
      x[i]=spin[rand()%spin.size()];
    }
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
    }*/
  int i = 0;
  cout << get<0>(matrix)[0] << endl;
  
  return 0;
}
