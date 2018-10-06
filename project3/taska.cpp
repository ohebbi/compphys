#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include<fstream>

using namespace std;

int main(int argc,char* argv[]){
    ofstream myfile;

    int n = 10000;

    vector<double> vx(n, 0);
    vector<double> vy(n, 0);
    vector<double> vz(n, 0);
    vector<double> rx(n, 0);
    vector<double> ry(n, 0);
    vector<double> rz(n, 0);
    vector<double> ax(n, 0);
    vector<double> ay(n, 0);

    vx[0]=(-1.24465);
    vy[0]=(6.132);
    vz[0]=(-0.00044);
    rx[0]=(1.0);
    ry[0]=(0.0);
    rz[0]=(0.0);
    
    float h = 1./n;

    if (strncmp(argv[1], "e", 2)==0){
      cout << "euler" << endl;
        for(int i = 1; i < n; i++){
            float r = sqrt(pow(rx[i-1], 2) + pow(ry[i-1], 2)+ pow(rz[i-1], 2));
           
            vx[i] = vx[i-1] - h * 4.0*pow(M_PI, 2.0) * rx[i-1] / pow(r,3.0);
            rx[i] = rx[i-1] + h * vx[i-1];

            vy[i] = vy[i-1] - h * 4.0*pow(M_PI, 2.0) * ry[i-1] / pow(r,3.0);
            ry[i] = ry[i-1] + h * vy[i-1];

            vz[i] = vz[i-1] - h * 4.0*pow(M_PI, 2.0) * rz[i-1] / pow(r,3.0);
            rz[i] = rz[i-1] + h * vz[i-1];

    }
    }
    else{
        cout << "Verlet bitches" << endl;
        for(int i = 1; i < n; i++){
            double r = sqrt(pow(rx[i-1],2)+pow(ry[i-1],2)+pow(rz[i-1],2));

            float ax = -4*pow(M_PI, 2.0)*rx[i-1]/pow(r, 3.0);

            rx[i] = rx[i-1]+h*vx[i-1]+pow(h,2.0)/2.0*ax;
            float ax1 = -4*pow(M_PI, 2.0)*rx[i]/pow(r, 3.0);
            vx[i] = vx[i-1]+h/2.0*(ax1+ax);


            float ay = -4*pow(M_PI, 2.0)*ry[i-1]/pow(r, 3.0);

            ry[i] = ry[i-1]+h*vy[i-1]+pow(h,2.0)/2.0*ay;
            float ay1 = -4*pow(M_PI, 2.0)*ry[i]/pow(r, 3.0);
            vy[i] = vy[i-1]+h/2.0*(ay1+ay);

            float az = -4*pow(M_PI, 2.0)*rz[i-1]/pow(r, 3.0);

            rz[i] = rz[i-1]+h*vz[i-1]+(pow(h,2.0)/2.0)*az;
            float az1 = -4*pow(M_PI, 2.0)*rz[i]/pow(r, 3.0);
            vz[i] = vz[i-1]+h/2.0*(az1+az);

    }
    cout << "Verlet out" << endl;
    
    }
    /*
    cout << "vx" <<" " << "vy" <<" " << "rx"<<" " << "ry" << endl;
    for(int i = 0; i< vy.size(); i++){
        cout << "vx: " << vx[i] << " vy: " << vy[i] << " rx: " << rx[i] << " ry " << ry[i] << endl;

    }

    
*/
    myfile.open("values3.txt");

    for (int i = 0; i < rx.size(); i++){
        myfile << rx[i] << " " << ry[i] << " " << rz[i] << endl;
    }
    myfile.close();
}
