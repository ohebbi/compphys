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


int bane(float final_time, Planet first, Planet second){
    int n = 1000000;

    double msun = 2e30;
    
    double vx = (first.inv[0]);
    double vy = (first.inv[1]);
    double vz = (first.inv[2]);
    
    double rx=(first.inv[3]);
    double ry=(first.inv[4]);
    double rz=(first.inv[5]);

    double jvx = (second.inv[0]);
    double jvy = (second.inv[1]);
    double jvz = (second.inv[2]);

    double jrx = (second.inv[3]);
    double jry = (second.inv[4]);
    double jrz = (second.inv[5]);
    
    double b = 3.0;
    double h = final_time/n;
    double r = sqrt(pow(rx,2)+pow(ry,2)+pow(rz,2));
    double ax1 = -4*pow(M_PI, 2.0)*rx/pow(r, b);
    double ay1 = -4*pow(M_PI, 2.0)*ry/pow(r, b);
    double az1 = -4*pow(M_PI, 2.0)*rz/pow(r, b);

    double jr = sqrt(pow(jrx, 2)+pow(jry,2)+pow(jrz,2));
    double jax1 = -4*pow(M_PI, 2.0)*jrx/pow(jr, b);
    double jay1 = -4*pow(M_PI, 2.0)*jry/pow(jr, b);
    double jaz1 = -4*pow(M_PI, 2.0)*jrz/pow(jr, b);

    double rej;

    ofstream tmpfile;
    tmpfile.open("values3.txt");
    
    for(int i = 1; i < n; i++){

        rej = pow((jrx-rx),2)+pow((jry-ry),2)+pow((jrz-rz),2);

	double ax = ax1;
        double ay = ay1;
        double az = az1;
     
        rx += h*vx+pow(h,2.0)/2.0*ax;
	ry += h*vy+pow(h,2.0)/2.0*ay;
	rz += h*vz+(pow(h,2.0)/2.0)*az;
	
	r = sqrt(pow(rx,2)+pow(ry,2)+pow(rz,2));
	
	ax1 = -4*pow(M_PI, 2.0)*(rx/pow(r, b)+first.mass/(rej*msun));
        vx += (h/2.0)*(ax1+ax);    
        
        ay1 = -4*pow(M_PI, 2.0)*(ry/pow(r, b)+first.mass/(rej*msun));
        vy += (h/2.0)*(ay1+ay);
       
        az1 = -4*pow(M_PI, 2.0)*(rz/pow(r, b)+first.mass/(rej*msun));
        vz += (h/2.0)*(az1+az);
	
	if( i%1000==0){
	  tmpfile << rx << " " << ry << " " << rz << "\n";
	}

	double jax = jax1;
        double jay = jay1;
        double jaz = jaz1;
     
        jrx += h*jvx+pow(h,2.0)/2.0*jax;
	jry += h*jvy+pow(h,2.0)/2.0*jay;
	jrz += h*jvz+(pow(h,2.0)/2.0)*jaz;
	
	jr = sqrt(pow(jrx,2)+pow(jry,2)+pow(jrz,2));
	
	jax1 = -4*pow(M_PI, 2.0)*(jrx/pow(jr, b)+second.mass/(rej*msun));
        jvx += (h/2.0)*(jax1+jax);    
        
        jay1 = -4*pow(M_PI, 2.0)*(jry/pow(jr, b)+second.mass/(rej*msun));
        jvy += (h/2.0)*(jay1+jay);
       
        jaz1 = -4*pow(M_PI, 2.0)*(jrz/pow(jr, b)+second.mass/(rej*msun));
        jvz += (h/2.0)*(jaz1+jaz);
	
	if( i%1000==0){
	  tmpfile << jrx << " " << jry << " " << jrz << "\n";
	}


}
    tmpfile.close();
    return 0;
}

int main(int argc,char* argv[]){
    
    double final_time = 20000.0;

    Planet earth;
    earth.mass = 6.0e24;
    earth.DS = 1;
    earth.inv = {-1.3497, 6.1340, 0.00042945, 0.980205, 0.205756, -0.0000874989};
    //earth.inv = {0, 2*M_PI, 0, 1, 0, 0};
  

    Planet jupiter;
    jupiter.mass = 1.9e27;
    jupiter.DS = 5.2;
    jupiter.inv = {2.342, -1.264, -0.04712, -2.71845, -4.6282558, 0.0800034};

    return bane(final_time, earth, jupiter);
}
