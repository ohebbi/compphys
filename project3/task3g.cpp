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

vector<double> force(Planet planets){
        double rej;
        vector<double> f = {0,0,0};
        double rej_2;
        double c2 = 9.0e16; //(speed of light)^2
       double r;
        double l;
        double v;
        
	r = sqrt(pow(planets.pos[0], 2)+pow(planets.pos[1],2)+pow(planets.pos[2],2));
        v = sqrt(pow(planets.pos[3], 2) + pow(planets.pos[4],2)+pow(planets.pos[5],2));
        rej = pow(r,3);
        l = pow((planets.pos[0]*planets.pos[4]-planets.pos[1]*planets.pos[3]),2); // l ^2 ; *sin \theta, but this is uno.
	//l = v*v;
        for(int j = 0; j<f.size(); j++){
	  f[j] = (planets.pos[j]/(rej))*(1+((3*l)/(r*r*c2)));
                }
        
        return f;
}


vector<double> instal(Planet first){      
        
        first.pos = first.inv;      
  
        vector<double> f = force(first);

        double r = sqrt(pow(first.pos[0],2)+pow(first.pos[1],2));
        double ax1 = -4*pow(M_PI, 2.0)*f[0];
        double ay1 = -4*pow(M_PI, 2.0)*f[1];
        double az1 = -4*pow(M_PI, 2.0)*f[2];

        return {first.pos[0], first.pos[1], first.pos[2], first.pos[3], first.pos[4], first.pos[5], ax1, ay1, az1};

}

vector<double> get_pos(Planet p, double h, double b){

        p.pos[0] += h*p.pos[3]+(pow(h,2.0)/2.0)*p.pos[6];
	p.pos[1] += h*p.pos[4]+(pow(h,2.0)/2.0)*p.pos[7];
	p.pos[2] += h*p.pos[5]+(pow(h,2.0)/2.0)*p.pos[8];

	vector<double> f = force(p);

	double r = sqrt(pow(p.pos[0],2)+pow(p.pos[1],2));

	double ax = -4*pow(M_PI, 2.0)*f[0];
        p.pos[3] += (h/2.0)*(ax+p.pos[6]);

        double ay = -4*pow(M_PI, 2.0)*f[1];
        p.pos[4] += (h/2.0)*(ay+p.pos[7]);

        double az = -4*pow(M_PI, 2.0)*f[2];
        p.pos[5] += (h/2.0)*(az+p.pos[8]);

        return {p.pos[0], p.pos[1], p.pos[2], p.pos[3], p.pos[4], p.pos[5], ax, ay, az};

}

int bane(float final_time, double b, Planet p){
  int n = 1e7;
    
    double h = final_time/n;

    p.pos = instal(p);
    
    ofstream tmpfile;
    tmpfile.open("values2.txt"); //ordinary plotting
    //tmpfile.open("values2.xyz"); //for ovito

    double theta;
    double tol = 1e-11;

  
    double farten;
    double dis;

    vector<double> peri = {1000, 1000};
    
    
    for(int ii = 1; ii < n; ii++){

      p.pos = get_pos(p, h, b);
      /*
        for(int t = 0; t < peri.size(); t++){
                if(peri[t] < p.pos[t]){
                        peri[t] = p.pos[t];
                }
        }
        if(ii % 100/0.2 == 0){
                tmpfile << atan(peri[1]/peri[0]) << " " << h*ii <<  "\n";
                peri = {1000, 1000};
        }
      */
	/*
        if(ii%1000==0){
          for (int jj = 0; jj < planets.size(); jj++){

                tmpfile <<  planets[jj].pos[0] << " " << planets[jj].pos[1] << " " << planets[jj].pos[2] << " " << jj << " " << "\n";
          }
 	}
        
	       
	*/
	//tmpfile << 180*(theta)/M_PI << " " << h*ii << " " << p.pos[1] << " " << p.pos[0] <<"\n";
      	
      farten = sqrt(pow(p.pos[3],2)+pow(p.pos[4],2)+pow(p.pos[5],2));
      dis = sqrt(pow(p.pos[0],2)+pow(p.pos[1],2)+pow(p.pos[2],2));
        if ((fabs(dis-0.3075000) <= tol) and (fabs(farten-12.44) <= tol)){
	        theta =  atan(p.pos[1]/p.pos[0]);
		tmpfile << (theta) << " " << h*ii << " " << p.pos[1] << " " << p.pos[0] <<"\n";
         
        }
	
        //if(ii%1000==0){

        //}


        

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

    double final_time = 100;
    double b = 3.0;

    Planet mercury;
    mercury.name = "mercury";
    mercury.mass = 3.3e23;
    mercury.inv = {0.30750000, 0,0,0,12.44000,0};

    return bane(final_time, b, mercury);
}
