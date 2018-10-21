#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include<fstream>
#include <time.h>

using namespace std;

int main(int argc,char* argv[]){
    ofstream myfile;

    int n = pow(10,3);

    vector<double> vx(n, 0);
    vector<double> vy(n, 0);
    vector<double> vz(n, 0);
    vector<double> rx(n, 0);
    vector<double> ry(n, 0);
    vector<double> rz(n, 0);
    vector<double> ax(n, 0);
    vector<double> ay(n, 0);

    //vx[0]=(-1.24465);
    //vy[0]=(6.132);
    //vz[0]=(-0.00044);

    //Initial value for a circular orbit
    vx[0]=(0);
    vy[0]=(2*M_PI);
    vz[0]=(0);
    rx[0]=(1.0);
    ry[0]=(0.0);
    rz[0]=(0.0);

    float h = 1./n;
    double start;
    if (strncmp(argv[1], "e", 2)==0){
      cout << "euler" << endl;
      clock_t start = clock();
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

	double r;
	double ax;
	double ay;
	double az;
	double ax1;
	double ay1;
	double az1;
	clock_t start = clock();
  double constant = -4*pow(M_PI, 2.0);
        for(int i = 1; i < n; i++){
            r = sqrt(pow(rx[i-1],2)+pow(ry[i-1],2)+pow(rz[i-1],2));

            ax = constant*rx[i-1]/pow(r, 3.0);
	          ay = constant*ry[i-1]/pow(r, 3.0);
	          az = constant*rz[i-1]/pow(r, 3.0);

            rx[i] = rx[i-1]+h*vx[i-1]+pow(h,2.0)/2.0*ax;
	          ry[i] = ry[i-1]+h*vy[i-1]+pow(h,2.0)/2.0*ay;
	          rz[i] = rz[i-1]+h*vz[i-1]+(pow(h,2.0)/2.0)*az;

	          r = sqrt(pow(rx[i],2)+pow(ry[i],2)+pow(rz[i],2));

            ax1 = -4*pow(M_PI, 2.0)*rx[i]/pow(r, 3.0);
            vx[i] = vx[i-1]+h/2.0*(ax1+ax);

            ay1 = -4*pow(M_PI, 2.0)*ry[i]/pow(r, 3.0);
            vy[i] = vy[i-1]+h/2.0*(ay1+ay);

            az1 = -4*pow(M_PI, 2.0)*rz[i]/pow(r, 3.0);
            vz[i] = vz[i-1]+h/2.0*(az1+az);

    }

    cout << "Verlet out" << endl;


    }
    clock_t end=clock();
    double time = (double) (end-start)/ CLOCKS_PER_SEC;
    cout << "Calculating time for n="<< n << ':'<< time << "s"<< '\n';
    cout << '\n' << endl;
    double start1 = pow(pow(rx[0],2)+pow(ry[0],2)+pow(rz[0],2), 0.5);
    double slutt = pow(pow(rx[n-1],2)+pow(ry[n-1],2)+pow(rz[n-1],2), 0.5);
    double svar = start1 - slutt;
    cout << "start = " << start1 << endl;
    cout << "slutt = " << slutt << endl;
    cout << "forskjell = " << svar << endl;
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
