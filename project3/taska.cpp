#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include<fstream>

using namespace std;

int main(int argc,char* argv[]){
    int n = 10000;

    ofstream myfile;
  
    vector<double> vx(n, 1);
    vector<double> vy(n, 1);
    vector<double> rx(n, 1);
    vector<double> ry(n, 1);

    vx[0]=(-1.136);
    vy[0]=(6.172);
    rx[0]=(1);
    ry[0]=(0.0);

    
    
    float h = 1./n;

    if ( strncmp(argv[1], "e", 2)==0){
      cout << "euler" << endl;
        for(int i = 1; i < n; i++){
	  float r = sqrt(pow(rx[i-1],2)+pow(ry[i-1],2));
	  //float r = 1.0;
            vx[i] = vx[i-1]-h*4.0*pow(M_PI, 2.0)*rx[i-1]/pow(r,3.0);
            rx[i] = rx[i-1]+h*vx[i-1];

            vy[i] = vy[i-1]-h*4.0*pow(M_PI, 2.0)*ry[i-1]/pow(r,3.0);
            ry[i] = ry[i-1]+h*vy[i-1];

    }
    }
    else{
        for(int i = 1; i < n; i++){
            double r = sqrt(pow(rx[i-1],2)+pow(ry[i-1],2));

            float ax = -4*pow(M_PI, 2.0)*rx[i-1]/pow(r, 3.0);
            
            rx[i] = rx[i-1]+h*vx[i-1]+pow(h,2.0)/2.0*ax; 
            float ax1 = -4*pow(M_PI, 2.0)*rx[i]/pow(r, 3.0);          
            vx[i] = vx[i-1]+h*ax+pow(h,2.0)/2.0*((ax1-ax)/h);
            
            
            float ay = -4*pow(M_PI, 2.0)*ry[i-1]/pow(r, 3.0);
            
            ry[i] = ry[i-1]+h*vy[i-1]+pow(h,2.0)/2.0*ay; 
            float ay1 = -4*pow(M_PI, 2.0)*ry[i]/pow(r, 3.0);          
            vy[i] = vy[i-1]+h*ay+pow(h,2.0)/2.0*((ay1-ay)/h);

    }
    }
    cout << "vx" <<" " << "vy" <<" " << "rx"<<" " << "ry" << endl;
    for(int i = 0; i< vy.size(); i++){
      cout << vx[i] << "     " << vy[i] <<"     " << rx[i]<<"     " << ry[i] << endl;
        
        
    }   
        
    myfile.open ("values3.txt");
    
    for(int i = 0; i<n; i+=1){
        myfile << rx[i] << " " << ry[i] << endl;
    }
    myfile.close();
   
    

    return 0;

}   
