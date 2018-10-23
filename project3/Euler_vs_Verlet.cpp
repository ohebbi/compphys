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
    double final_time = 1.;
    int n = 1e7;

    //creating vectors to keep all information of the position, velocity and acceleration for earth
    vector<double> vx(n, 0);
    vector<double> vy(n, 0);
    vector<double> vz(n, 0);
    vector<double> rx(n, 0);
    vector<double> ry(n, 0);
    vector<double> rz(n, 0);
    vector<double> ax(n, 0);
    vector<double> ay(n, 0);
    double r;

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

    double h = final_time/n;
    double start;

    double hfopow = h * 4.0*pow(M_PI, 2.0);

    //to test the euler vs verlet, you can either run the program with "euler" in the command line, to run the euler method, or something else to run verlet
    if (strncmp(argv[1], "euler", 2)==0){ 
      
        cout << "euler" << endl;
        clock_t start = clock();//this starts the clock if you chose the euler method.
        for(int i = 1; i < n; i++){
            r = sqrt(pow(rx[i-1], 2) + pow(ry[i-1], 2)+ pow(rz[i-1], 2));

            vx[i] = vx[i-1] - hfopow * rx[i-1] / pow(r,3.0);
            rx[i] = rx[i-1] + h * vx[i-1];

            vy[i] = vy[i-1] - hfopow * ry[i-1] / pow(r,3.0);
            ry[i] = ry[i-1] + h * vy[i-1];

            vz[i] = vz[i-1] - hfopow * rz[i-1] / pow(r,3.0);
            rz[i] = rz[i-1] + h * vz[i-1];

    }
    }
    else{
        cout << "Verlet" << endl;

	
	    double ax;
	    double ay;
	    double az;
	    double ax1;
	    double ay1;
	    double az1;
        double h2 = h/2.0;
	    clock_t start = clock(); //this starts the clock if you chose the verlet method
    
        double fopow = -4*pow(M_PI, 2.0);
            for(int i = 1; i < n; i++){
            r = sqrt(pow(rx[i-1],2)+pow(ry[i-1],2)+pow(rz[i-1],2));

            ax = fopow*rx[i-1]/pow(r, 3.0);
	        ay = fopow*ry[i-1]/pow(r, 3.0);
	        az = fopow*rz[i-1]/pow(r, 3.0);

            rx[i] = rx[i-1]+h*vx[i-1]+pow(h,2.0)/2.0*ax;
	        ry[i] = ry[i-1]+h*vy[i-1]+pow(h,2.0)/2.0*ay;
	        rz[i] = rz[i-1]+h*vz[i-1]+(pow(h,2.0)/2.0)*az;

	        r = sqrt(pow(rx[i],2)+pow(ry[i],2)+pow(rz[i],2));

            ax1 = fopow*rx[i]/pow(r, 3.0);
            vx[i] = vx[i-1]+h*(ax1+ax);

            ay1 = fopow*ry[i]/pow(r, 3.0);
            vy[i] = vy[i-1]+h2*(ay1+ay);

            az1 = fopow*rz[i]/pow(r, 3.0);
            vz[i] = vz[i-1]+h2*(az1+az);

        }

        cout << "Verlet out" << endl;


    }
    clock_t end=clock();
    double time = (double) (end-start)/CLOCKS_PER_SEC;
    cout << "Calculating time for n="<< n << ':'<< time << "s"<< '\n';
    cout << '\n' << endl;
    double start1 = pow(pow(rx[0],2)+pow(ry[0],2)+pow(rz[0],2), 0.5);
    double slutt = pow(pow(rx[n-1],2)+pow(ry[n-1],2)+pow(rz[n-1],2), 0.5);
    double svar = start1 - slutt;
    cout << "start = " << start1 << endl;
    cout << "slutt = " << slutt << endl;
    cout << "forskjell = " << svar << endl;

    myfile.open("values3.txt");

    for (int i = 0; i < rx.size(); i++){
        if(i % 10000 == 0){
            myfile << rx[i] << " " << ry[i] << " " << rz[i] << " \n";
        }
        
    }
    myfile.close();
}
