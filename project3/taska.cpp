#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdio.h>

using namespace std;

int main(int argc,char* argv[]){
    int n = 100;

    vector<double> vx(n, 0);
    vector<double> vy(n, 0);
    vector<double> rx(n, 0);
    vector<double> ry(n, 0);
    vector<double> ax(n, 0);
    vector<double> ay(n, 0);

    vx[0]=(-1/sqrt(2.0*3.14159));
    vy[0]=(1/sqrt(2.0*3.14159));
    rx[0]=(1.0);
    ry[0]=(0.0);
    ax[0]=(4.0 * pow(3.14159, 2.0 )
    ay[0]=



    float h = 1./n;

    if ( strncmp(argv[1], "e", 2)==0){
      cout << "euler" << endl;
        for(int i = 1; i < n; i++){
            //float r = sqrt(pow(rx[i-1], 2) + pow(ry[i-1], 2));
            float r = 1.0;
            vx[i] = vx[i-1] + h * 4.0*pow(3.14159, 2.0) * rx[i-1] / pow(r,3.0);
            rx[i] = rx[i-1] + h * vx[i-1];

            vy[i] = vy[i-1] + h * 4.0*pow(3.14159, 2.0) * ry[i-1] / pow(r,3.0);
            ry[i] = ry[i-1] + h * vy[i-1];

    }
    }
    else{
        for(int i = 1; i < n; i++){
            double r = sqrt(pow(rx[i],2)+pow(ry[i],2));

            float ax = -4*pow(M_PI, 2.0)*rx[i-1]/pow(r, 3.0);

            rx[i] = rx[i-1]+h*vx[i]+pow(h,2.0)/2.0*ax;
            float ax1 = -4*pow(M_PI, 2.0)*rx[i]/pow(r, 3.0);
            vx[i] = vx[i-1]+h*ax+pow(h,2.0)/2.0*((ax1-ax)/h);


            float ay = -4*pow(M_PI, 2.0)*ry[i-1]/pow(r, 3.0);

            ry[i] = ry[i-1]+h*vy[i]+pow(h,2.0)/2.0*ay;
            float ay1 = -4*pow(M_PI, 2.0)*ry[i]/pow(r, 3.0);
            vy[i] = vy[i-1]+h*ay+pow(h,2.0)/2.0*((ay1-ay)/h);

    }
    }
    for(int i = 0; i< vy.size(); i++){
        cout << "vx: " << vx[i] << " vy: " << vy[i] << " rx: " << rx[i] << " ry " << ry[i] << endl;

    }



}
