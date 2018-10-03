#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

int main(int argc,char* argv[]){

    vector<double> vx;
    vector<double> vy;
    vector<double> rx;
    vector<double> ry;

    vx.push_back(1/sqrt(2));
    vy.push_back(1/sqrt(2));
    rx.push_back(1);
    ry.push_back(0);

    
    int n = 100;
    float h = 1./n;

    if (argv[1]= "euler"){
        for(int i = 1; i < n; i++){
            double r = sqrt(pow(rx[i],2)+pow(ry[i],2));

            vx[i] = vx[i-1]-h*4*pow(pi, 2.0)*rx[i-1]/pow(r,3);
            rx[i] = rx[i-1]+h*vx[i];

            vy[i] = vy[i-1]-h*4*pow(pi, 2.0)*ry[i-1]/pow(r,3);
            ry[i] = ry[i-1]+h*vy[i];

    }
    }
    else{
        for(int i = 1; i < n; i++){
            double r = sqrt(pow(rx[i],2)+pow(ry[i],2));

            float ax = -4*pow(pi, 2.0)*rx[i-1]/pow(r, 3.0);
            
            rx[i] = rx[i-1]+h*vx[i]+pow(h,2.0)/2.0*ax; 
            float ax1 = -4*pow(pi, 2.0)*rx[i]/pow(r, 3.0);          
            vx[i] = vx[i-1]+h*ax+pow(h,2.0)/2.0*((ax1-ax)/h);
            
            
            float ay = -4*pow(pi, 2.0)*ry[i-1]/pow(r, 3.0);
            
            ry[i] = ry[i-1]+h*vy[i]+pow(h,2.0)/2.0*ay; 
            float ay1 = -4*pow(pi, 2.0)*ry[i]/pow(r, 3.0);          
            vy[i] = vy[i-1]+h*ay+pow(h,2.0)/2.0*((ay1-ay)/h);

    }

    for(int i = 0; i< vy.size(); i++){
        cout << "vx: " << vx[i] << endl;
        cout << "vy: " << vy[i] << endl;
        cout << "rx: " << rx[i] << endl;
        cout << "ry " << ry[i] << endl;
        
    }   
        
    }

}   
