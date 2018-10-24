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
    myfile.open("values3.txt");

    double final_time = 10.; //gives the number of years the simulation shall run over
    int n = 1e9; //number of time steps
    double h = final_time/n; //time step

    //creating vectors to keep all information of the position, velocity and acceleration for earth
    vector<double> vx(2, 0);
    vector<double> vy(2, 0);
    vector<double> vz(2, 0);
    vector<double> rx(2, 0);
    vector<double> ry(2, 0);
    vector<double> rz(2, 0);
    vector<double> ax(2, 0);
    vector<double> ay(2, 0);
    vector<double> az(2, 0);
    double r;

    //initial values
    rx[0]=(1.0);
    ry[0]=(0.0);
    rz[0]=(0.0);
    //vx[0]=(-1.24465);
    //vy[0]=(6.132);
    //vz[0]=(-0.00044);

    //Initial value for a circular orbit
    vx[0]=(0);
    vy[0]=(2*M_PI);
    vz[0]=(0);

    double start1 = sqrt(pow(rx[0],2)+pow(ry[0],2)+pow(rz[0],2));
    double slutt;
    
    double start;

    //creating constants outside the loop to reduce the number of flops
    double hfopow = h * 4.0*pow(M_PI, 2.0);
    double rhfopow;

    //to test the euler vs verlet, you can either run the program with "euler" in the command line, to run the euler method, or something else to run verlet
    if (strncmp(argv[1], "euler", 2)==0){

        cout << "euler" << endl;
        clock_t start = clock();//this starts the clock if you chose the euler method.

        for(int i = 1; i < n; i++){//runs the Euler method

            rhfopow = hfopow/pow(sqrt(pow(rx[0], 2)+ pow(ry[0], 2)+ pow(rz[0], 2)),3);

            vx[1] = vx[0]-rhfopow*rx[0];
            rx[1] = rx[0]+h*vx[0];
            vx[0] = vx[1];
            rx[0] = rx[1];

            vy[1] = vy[0]-rhfopow*ry[0];
            ry[1] = ry[0]+h*vy[0];
            vy[0] = vy[1];
            ry[0] = ry[1];

            vz[1] = vz[0]-rhfopow*rz[0];
            rz[1] = rz[0]+h*vz[0];
            vz[0] = vz[1];
            rz[0] = rz[1];

            if(i % 10000 == 0){
                myfile << rx[0] << " " << ry[0] << " " << rz[0] << " \n";//coordinates shall be read from a python program
            }

        }
	slutt = sqrt(pow(rx[0],2)+pow(ry[0],2)+pow(rz[0],2));
	myfile.close();
    }


    else{
        cout << "Verlet" << endl;

        //initial values installed from the vectors. We don't want to use vectors for position and velocity in this case.
        double rxv = rx[0];
        double ryv = ry[0];
        double rzv = rz[0];
        double vxv =vx[0];
        double vyv = vy[0];
        double vzv = vz[0];

        //defining constants outside the loop to reduce the number of flops
        double h2 = h/2.0;
        double hpow2 = pow(h,2.0)/2.0;
        r = pow(sqrt(pow(rxv,2)+pow(ryv,2)+pow(rzv,2)),3);
        double fopow = -4*pow(M_PI, 2.0);
        double rfopow = fopow/r;

        ax[0] = rfopow*rxv;
	    ay[0] = rfopow*ryv;
	    az[0] = rfopow*rzv;
	    clock_t start = clock(); //this starts the clock if you chose the verlet method


        for(int i = 1; i < n; i++){//runs the Verlet method

	  rxv += h*vxv+hpow2*ax[0];
	        ryv += h*vyv+hpow2*ay[0];
	        rzv += h*vzv+hpow2*az[0];

	        rfopow = fopow/pow(sqrt(pow(rxv,2)+pow(ryv,2)+pow(rzv,2)),3);

            ax[1] = rfopow*rxv;
            vxv += h2*(ax[1]+ax[0]);
            ax[0] = ax[1];

            ay[1] = rfopow*ryv;
            vyv += h2*(ay[1]+ay[0]);
            ay[0] = ay[1];

            az[1] = rfopow*rzv;
            vzv += h2*(az[1]+az[0]);
            az[0] = az[1];

            if(i % 10000 == 0){
                myfile << rxv << " " << ryv << " " << rzv << " \n";//coordinates shall be read from a python program
            }
        }
	slutt = sqrt(pow(rxv,2)+pow(ryv,2)+pow(rzv,2));
        myfile.close();
    }
    clock_t end=clock();

    //prints the time used
    double time = (double) (end-start)/CLOCKS_PER_SEC;
    cout << "Calculating time for n="<< n << ':'<< time << "s"<< '\n';
    cout << '\n' << endl;
    
    
    double svar = start1 - slutt;
    cout << "start = " << start1 << endl;
    cout << "slutt = " << slutt << endl;
    cout << "forskjell = " << svar << endl;
    

    




}
