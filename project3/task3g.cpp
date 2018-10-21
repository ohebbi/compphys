#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include<fstream>

using namespace std;

class Planet {
    public:
        string name;
        double mass;
        double DS;
        vector<double> inv;
        vector<double> pos;

}globalPlanet;

vector<double> force(vector<Planet> planets, double msun){
        double rej;
        vector<double> f = {0,0,0};
        double rej_2;
        double c2 = 9.0e16; //(speed of light)^2
        double r;
        double l;
        double v;

        for(int i = 1; i < planets.size(); i++){
                r = sqrt( pow((planets[i].pos[0] - planets[0].pos[0]), 2) + pow( (planets[i].pos[1] - planets[0].pos[1]),2) + pow((planets[i].pos[2] - planets[0].pos[2]), 2));
                v = sqrt( pow((planets[i].pos[3] - planets[0].pos[3]), 2) + pow( (planets[i].pos[4] - planets[0].pos[4]),2) + pow((planets[i].pos[5] - planets[0].pos[5]), 2));
                rej = pow(r,3);
                l = r*r*v*v; // l ^2 ; *sin \theta, but this is uno.

                for(int j = 0; j<f.size(); j++){
                        f[j] += ((planets[0].pos[j]-planets[i].pos[j])*planets[i].mass/(rej*msun))*(1+((3*l)/(r*r*c2)));

                }
        }
        return f;
}


vector<double> instal(double b, vector<Planet> planets, int in){

        swap(planets[0], planets[in]);


        double msun = 1.989e30;
        Planet first = planets[0];
        for (int i = 0; i< planets.size(); i++){
                planets[i].pos = planets[i].inv;
        }

        vector<double> f = force(planets, msun);

        double r = sqrt(pow(first.inv[3],2)+pow(first.inv[4],2)+pow(first.inv[5],2));
        double ax1 = -4*pow(M_PI, 2.0)*(first.inv[3]/pow(r, b)+f[0]);
        double ay1 = -4*pow(M_PI, 2.0)*(first.inv[4]/pow(r, b)+f[1]);
        double az1 = -4*pow(M_PI, 2.0)*(first.inv[5]/pow(r, b)+f[2]);

        return {first.inv[3], first.inv[4], first.inv[5], first.inv[0], first.inv[1], first.inv[2], ax1, ay1, az1};

}

vector<double> get_pos(vector<Planet> planets, double h, double b, int in){
        swap(planets[0], planets[in]);

        Planet p = planets[0];
        double msun = 2e30;

        p.pos[0] += h*p.pos[3]+(pow(h,2.0)/2.0)*p.pos[6];
	      p.pos[1] += h*p.pos[4]+(pow(h,2.0)/2.0)*p.pos[7];
	      p.pos[2] += h*p.pos[5]+(pow(h,2.0)/2.0)*p.pos[8];

	      vector<double> f = force(planets, msun);

	      double r = sqrt(pow(p.pos[0],2)+pow(p.pos[1],2)+pow(p.pos[2],2));

	      double ax = -4*pow(M_PI, 2.0)*(p.pos[0]/pow(r, b)+f[0]);
        p.pos[3] += (h/2.0)*(ax+p.pos[6]);

        double ay = -4*pow(M_PI, 2.0)*(p.pos[1]/pow(r, b)+f[1]);
        p.pos[4] += (h/2.0)*(ay+p.pos[7]);

        double az = -4*pow(M_PI, 2.0)*(p.pos[2]/pow(r, b)+f[2]);
        p.pos[5] += (h/2.0)*(az+p.pos[8]);

        return {p.pos[0], p.pos[1], p.pos[2], p.pos[3], p.pos[4], p.pos[5], ax, ay, az};

}

int bane(float final_time, double b, vector<Planet> planets){
    int n = 10e6;

    double h = final_time/n;

    planets[0].pos = instal(b, planets, 0);




    ofstream tmpfile;
    tmpfile.open("values2.txt"); //ordinary plotting
    //tmpfile.open("values2.xyz"); //for ovito

    double theta;
    double tol = 10e-5;

    double lollipopp;
    double farten;


    for(int ii = 1; ii < n; ii++){

        planets[0].pos = get_pos(planets, h, b, 0);


        farten = pow(pow(planets[0].pos[3]+planets[0].pos[4]+planets[0].pos[5],2),0.5);
        if (planets[0].pos[0]>=0 and (planets[0].pos[1]>=0)){
          if (farten < (12.44+tol)){
            if (farten > (12.44-tol)){
              theta = atan(planets[0].pos[1]/planets[0].pos[0]);
              cout << farten << " " << theta*3600 << " " << h*ii << " " << endl;

              tmpfile << planets[0].pos[0] << " " << planets[0].pos[1] << " " << planets[0].pos[2] << " " << "\n";

            }
          }
        //if(ii%1000==0){

        //}


        }

          //for ordinary plotting; comment out the two next lines
          //tmpfile << "1" << "\n"; //number of planets
          //tmpfile << "commentline that needs to be here" << "balle" << "\n";



                //plotting in ovito
                //tmpfile << jj << " " << planets[jj].pos[0] << " " << planets[jj].pos[1] << " " << planets[jj].pos[2] << " " << jj << " " << "\n";








}
    tmpfile.close();
    return 0;
}

int main(int argc,char* argv[]){

    double final_time = 150;
    double b = 3.0;

    Planet sun;
    sun.name = "sun";
    sun.mass = 1.989e30;
    sun.DS = 0;
    sun.inv={0,0,0, 1,0,0};

    Planet earth;
    earth.name = "earth";
    earth.mass = 6.0e24;
    earth.DS = 1;
    earth.inv = {-1.3497, 6.1340, 0.00042945, 0.980205, 0.205756, -0.0000874989};
    //earth.inv = {0, 2*M_PI, 0, 1, 0, 0};

    Planet jupiter;
    jupiter.name = "jupiter";
    jupiter.mass = 1.9e27;
    jupiter.DS = 5.2;
    jupiter.inv = {2.342, -1.264, -0.04712, -2.71845, -4.6282558, 0.0800034};

    Planet mars;
    mars.name = "mars";
    mars.mass = 6.6e23;
    mars.DS = 1.52;
    mars.inv = {1.318, 5.420, 0.081, 1.3490045, -0.29756, -0.03956};

    Planet mercury;
    mercury.name = "mercury";
    mercury.mass = 3.3e23;
    //mercury.DS = 0.39;
    //mercury.inv = {4.82, -7.13, -1.025, -0.3294, -0.29288, 0.005618};
    mercury.DS = 0.3075;
    //mercury.inv = {0.3075,0,0,0,12.44,0};
    mercury.inv = {0,12.44,0,0.3075,0,0};


    Planet venus;
    venus.name = "venus";
    venus.mass = 4.9e24;
    venus.DS = 0.72;
    venus.inv = {0.37, 7.34, 0.079, 0.7244, -0.03279, -0.042427};

    Planet saturn;
    saturn.name = "saturn";
    saturn.mass = 5.5e26;
    saturn.DS = 9.54;
    saturn.inv = {1.901, 0.296, -0.0807, 1.497, -9.9435, 0.1133};

    Planet uranus;
    uranus.name = "uranus";
    uranus.mass = 8.8e25;
    uranus.DS = 19.19;
    uranus.inv = {-0.73, 1.175, 0.0139, 17.196, 0.9965, -0.1858};

    Planet neptune;
    neptune.name = "neptune";
    neptune.mass = 1.03e26;
    neptune.DS = 30.06;
    neptune.inv = {0.289, 1.11, -0.0297, 28.91, -7.753, -0.5067};

    Planet pluto;
    pluto.name = "pluto";
    pluto.mass = 1.31e22;
    pluto.DS = 39.53;
    pluto.inv = {1.0975, 0.1535, -0.331, 11.614, -31.58, 0.01979};

    vector<Planet> planets =  {mercury};

    return bane(final_time, b, planets);
}
