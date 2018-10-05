#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include<fstream>



using namespace std;


class Planet {
    public:
        double mass;
        double DS;
        vector<double> inv;
        
        
}globalPlanet;



int bane(float final_time, vector<double> init){
    int n = 10000000;

    double vx = (init[0]);
    double vy = (init[1]);
    double vz = (init[2]);
    
    double rx=(init[3]);
    double ry=(init[4]);
    double rz=(init[5]);
    
    double b = 3.;
    double h = final_time/n;
    double r;
    r = sqrt(pow(rx,2)+pow(ry,2)+pow(rz,2));
    double ax1 = -4*pow(M_PI, 2.0)*rx/pow(r, b);
    double ay1 = -4*pow(M_PI, 2.0)*ry/pow(r, b);
    double az1 = -4*pow(M_PI, 2.0)*rz/pow(r, b);

    ofstream tmpfile;
    tmpfile.open("values3.txt");
    
    for(int i = 1; i < n; i++){
        
      
        double ax = ax1;
        double ay = ay1;
        double az = az1;
     

        rx += h*vx+pow(h,2.0)/2.0*ax;
	ry += h*vy+pow(h,2.0)/2.0*ay;
	rz += h*vz+(pow(h,2.0)/2.0)*az;
	
	r = sqrt(pow(rx,2)+pow(ry,2)+pow(rz,2));
	
	ax1 = -4*pow(M_PI, 2.0)*rx/pow(r, b);
        vx += (h/2.0)*(ax1+ax);
       
        
        ay1 = -4*pow(M_PI, 2.0)*ry/pow(r, b);
        vy += (h/2.0)*(ay1+ay);

        
        az1 = -4*pow(M_PI, 2.0)*rz/pow(r, b);
        vz += (h/2.0)*(az1+az);
	
	if( i%1000==0){
	  tmpfile << rx << " " << ry << " " << rz << "\n";
	}

}
    tmpfile.close();
    return 0;
}

int main(int argc,char* argv[]){
    
    double final_time = 1.0;

    Planet earth;
    earth.mass = 1;
    earth.DS = 1;
    earth.inv = {-1.3497, 6.1340, 0.00042945, 0.980205, 0.205756, -0.0000874989};
    //earth.inv = {0, 2*M_PI, 0, 1, 0, 0};
  

    //Planet jupiter;
    //jupiter.mass = 6.0e3/1.9;
    //jupiter.DS = 5.2;
    //jupiter.inv = {2.342, -1.264, -0.04712, -2.71845, -4.6282558, 0.0800034};

    return bane(final_time, earth.inv);
}
