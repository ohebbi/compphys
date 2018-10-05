#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include<fstream>
#include<tuple>


using namespace std;


class Planet {
    public:
        double mass;
        double DS;
        vector<double> inv;
        tuple<vector<double>, vector<double>, vector<double>> orbit;
        
}globalPlanet;



tuple<vector<double>, vector<double>, vector<double>> bane(float final_time, vector<double> init){
    int n = 100000;
    
    vector<double> vx(n, 0);
    vector<double> vy(n, 0);
    vector<double> vz(n, 0);
    vector<double> rx(n, 0);
    vector<double> ry(n, 0);
    vector<double> rz(n, 0);
    //vector<double> ax(n, 0);
    //vector<double> ay(n, 0);

    vx[0]=(init[0]);
    vy[0]=(init[1]);
    vz[0]=(init[2]);
    rx[0]=(init[3]);
    ry[0]=(init[4]);
    rz[0]=(init[5]);
    
    double b = 3.;
    double h = final_time/n;
    double r;
    r = sqrt(pow(rx[0],2)+pow(ry[0],2)+pow(rz[0],2));
    double ax1 = -4*pow(M_PI, 2.0)*rx[0]/pow(r, b);
    double ay1 = -4*pow(M_PI, 2.0)*ry[0]/pow(r, b);
    double az1 = -4*pow(M_PI, 2.0)*rz[0]/pow(r, b);
    
    
    for(int i = 1; i < n; i++){
        
      
        double ax = ax1;
        double ay = ay1;
        double az = az1;
     

        rx[i] = rx[i-1]+h*vx[i-1]+pow(h,2.0)/2.0*ax;
	ry[i] = ry[i-1]+h*vy[i-1]+pow(h,2.0)/2.0*ay;
	rz[i] = rz[i-1]+h*vz[i-1]+(pow(h,2.0)/2.0)*az;
	
	r = sqrt(pow(rx[i],2)+pow(ry[i],2)+pow(rz[i],2));
	
	ax1 = -4*pow(M_PI, 2.0)*rx[i]/pow(r, b);
        vx[i] = vx[i-1]+(h/2.0)*(ax1+ax);
       
        
        ay1 = -4*pow(M_PI, 2.0)*ry[i]/pow(r, b);
        vy[i] = vy[i-1]+(h/2.0)*(ay1+ay);

        
        az1 = -4*pow(M_PI, 2.0)*rz[i]/pow(r, b);
        vz[i] = vz[i-1]+(h/2.0)*(az1+az);


}
    return {rx, ry, rz};
}

int main(int argc,char* argv[]){
    ofstream myfile;
    double final_time = 1.0

    Planet earth;
    earth.mass = 1;
    earth.DS = 1;
    earth.inv = {-1.3497, 6.1340, 0.00042945 , 0.980205, 0.205756, -0.0000874989};
    //earth.inv = {0, 2*M_PI, 0 , 1, 0, 0};
    auto over = bane(final_time, earth.inv);
    earth.orbit = over;

    Planet jupiter;
    jupiter.mass = 6.0e3/1.9;
    jupiter.DS = 5.2;
    jupiter.inv = {2.342, -1.264, -0.04712, -2.71845, -4.6282558, 0.0800034};
    auto j_over = bane(final_time, jupiter.inv);


    myfile.open("values3.txt");

    for (int i = 0; i < get<0>(earth.orbit).size(); i++){
        myfile << get<0>(earth.orbit)[i] << " " << get<1>(earth.orbit)[i] << " " << get<2>(earth.orbit)[i] << "\n";
    }
    myfile.close();

}
