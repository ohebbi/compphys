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

vector<double> find_center(vector<Planet> planets){

        vector<double> center {0,0,0,0,0,0};
        double tot_mass = 0;

        for (int i = 0; i < planets.size(); i++){

                for (int j = 3; j < center.size(); j++){
                        center[j] -= planets[i].mass*planets[i].inv[j];
                }
                for (int jj = 0; jj < 3; jj++){
                        center[jj] -= planets[i].mass*planets[i].inv[jj]/planets[0].mass;
                }

                tot_mass += planets[i].mass;
        }

        for (int ii = 0; ii < center.size();ii++){
                center[ii] /= tot_mass;
        }
        return center;
}

vector<double> force(vector<Planet> planets, double msun){
        double rej;
        vector<double> f = {0,0,0};
        double rej_2;
        double c2 = 9.0e16; //(speed of light)^2
        double r;
        double l;
        double v;

        for(int i = 1; i < planets.size(); i++){
                r = sqrt(pow((planets[i].pos[0]-planets[0].pos[0]),2)+pow((planets[i].pos[1]-planets[0].pos[1]),2)+pow((planets[i].pos[2]-planets[0].pos[2]),2));
                v = sqrt(pow((planets[i].pos[3]-planets[0].pos[3]),2)+pow((planets[i].pos[4]-planets[0].pos[4]),2)+pow((planets[i].pos[5]-planets[0].pos[5]),2));
                rej = pow(r,3);
                l = r*r*v*v; //*sin \theta, but this is uno.

                for(int j = 0; j<f.size(); j++){
                        f[j] += ((planets[0].pos[j]-planets[i].pos[j])*planets[i].mass/(rej*msun))*(1+(3*l)/(r*r*c2));

                }
        }
        return f;
}

vector<double> instal(double b, vector<Planet> planets, int in, double msun){

        swap(planets[0], planets[in]);

        Planet first = planets[0];
        for (int i = 0; i< planets.size(); i++){
                planets[i].pos = planets[i].inv;
        }

        vector<double> f = force(planets, msun);

        double ax1 = -4*pow(M_PI, 2.0)*(f[0]);
        double ay1 = -4*pow(M_PI, 2.0)*(f[1]);
        double az1 = -4*pow(M_PI, 2.0)*(f[2]);

        return {first.inv[3], first.inv[4], first.inv[5], first.inv[0], first.inv[1], first.inv[2], ax1, ay1, az1};
}

vector<double> get_pos(vector<Planet> planets, double h, double b, int in, double msun){
        swap(planets[0], planets[in]);

        Planet p = planets[0];

        p.pos[0] += h*p.pos[3]+(pow(h,2.0)/2.0)*p.pos[6];
	      p.pos[1] += h*p.pos[4]+(pow(h,2.0)/2.0)*p.pos[7];
	      p.pos[2] += h*p.pos[5]+(pow(h,2.0)/2.0)*p.pos[8];

        vector<double> f = force(planets, msun);

	      double ax = -4*pow(M_PI, 2.0)*f[0];
        p.pos[3] += (h/2.0)*(ax+p.pos[6]);

        double ay = -4*pow(M_PI, 2.0)*f[1];
        p.pos[4] += (h/2.0)*(ay+p.pos[7]);

        double az = -4*pow(M_PI, 2.0)*(f[2]);
        p.pos[5] += (h/2.0)*(az+p.pos[8]);

        return {p.pos[0], p.pos[1], p.pos[2], p.pos[3], p.pos[4], p.pos[5], ax, ay, az};
}

int bane(float final_time, double b, vector<Planet> planets){
    int n = 1000000;
    double h = final_time/n;

    for(int i = 0; i < planets.size(); i++){
                planets[i].pos = instal(b, planets, i, 2e30);

        }

    ofstream tmpfile;
    tmpfile.open("values3.txt");

    for(int ii = 0; ii < n; ii++){


        for (int j = 0; j < planets.size(); j++){
                 planets[j].pos = get_pos(planets, h, b, j, 2e30);
        }

        if(ii%100==0){
          for (int jj = 0; jj < planets.size(); jj++){

                tmpfile <<  planets[jj].pos[0] << " " << planets[jj].pos[1] << " " << planets[jj].pos[2] << " " << jj << " " << "\n";
          }
 	}
}
    tmpfile.close();
	return 0;
}
/*    
    ofstream tmpfile;
    tmpfile.open("values2.txt");



    for(int ii = 1; ii < n; ii++){

        for (int j = 0; j < planets.size(); j++){
                 planets[j].pos = get_pos(planets, h, b, j,2e30);

        }

        if(ii%100==0){
          //for ordinary plotting; comment out the two next lines
          tmpfile << "2" << "\n"; //number of planets
          tmpfile << "commentline that needs to be here" << "balle" << "\n";

          for (int jj = 0; jj < planets.size(); jj++){
                tmpfile << jj <<" " << planets[jj].pos[0] << " " << planets[jj].pos[1] << " " << planets[jj].pos[2] << " " << jj << " " << "\n";
          }
  }
}
    tmpfile.close();
    return 0;
}
*/
int main(int argc,char* argv[]){
    double final_time = 100;
    double b = 3.0;

    Planet mercury;
    mercury.name = "mercury";
    mercury.mass = 3.3e23;
    mercury.DS = 0.39;
    mercury.inv = {4.82, -7.13, -1.025, -0.3294, -0.29288, 0.005618};

    Planet sola;
    sola.name = "soool";
    sola.mass = 2e30;
    sola.inv = {0,0,0,0,0,0};
    vector<Planet> planets = {sola, mercury};
    sola.inv = find_center(planets);

    planets[0] = sola;


    return bane(final_time, b, planets);
}
