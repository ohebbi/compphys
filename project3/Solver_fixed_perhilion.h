#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include<fstream>

using namespace std;


//here the decomposed force on a planet is calculated
vector<double> force(Planet* planets, double msun){
       double rej;
        vector<double> f = {0,0,0};
        double rej_2;
        double c2 = pow(63284.9,2); //(speed of light)^2
        double r;
        double l;
        double v;
        
	    r = sqrt(pow(planets->pos[0], 2)+pow(planets->pos[1],2)+pow(planets->pos[2],2));
        v = sqrt(pow(planets->pos[3], 2) + pow(planets->pos[4],2)+pow(planets->pos[5],2));
        rej = pow(r,3);
	    l = pow(0.3075*12.44,2); // l ^2 ; *sin \theta, but this is uno.
	
        for(int j = 0; j<f.size(); j++){
	        f[j] = (planets->pos[j]/(rej))*(1+((3*l)/(r*r*c2)));
            }
        
        return f;
}

//instal creates the initial acceleration needed for the Verlet method
vector<double> instal(vector<double> C, Planet* planets){
    
        
        Planet* first = planets;

        planets->pos = planets->inv;
        
        vector<double> f = force(planets, C[5]);
        double r = pow(sqrt(pow(first->pos[0],2)+pow(first->pos[1],2)+pow(first->pos[2],2)),C[0]);

        double ax1 = C[2]*f[0];
        double ay1 = C[2]*f[1];
        double az1 = C[2]*f[2];

        return {first->pos[0], first->pos[1], first->pos[2], first->pos[3], first->pos[4], first->pos[5], ax1, ay1, az1};
}

//this function returns both position, velocity and acceleration for the given planet
vector<double> get_pos(Planet* planets, vector<double> C){
        Planet* p = planets;
        
        for(int i = 0; i < 3; i++){
                p->pos[i] += C[1]*p->pos[i+3]+C[3]*p->pos[i+6];
        }
        
        vector<double> f = force(planets, C[5]);
        double r = pow(sqrt(pow(p->pos[0],2)+pow(p->pos[1],2)+pow(p->pos[2],2)),C[0]);

	    double ax = C[2]*f[0];
        p->pos[3] += C[4]*(ax+p->pos[6]);

        double ay = C[2]*f[1];
        p->pos[4] += C[4]*(ay+p->pos[7]);

        double az = C[2]*f[2];
        p->pos[5] += C[4]*(az+p->pos[8]);

        return {p->pos[0], p->pos[1], p->pos[2], p->pos[3], p->pos[4], p->pos[5], ax, ay, az};
}

//this the class that will solve our differential equation
class Solver_perihilion{
    
    public:
            int solution(float final_time, double b, vector<Planet*> planets, int n){
                
                Planet* p = planets[0]; 
                double h = final_time/n;     
                double hpo = pow(h,2.0)/2.0;
                double hdiv = (h/2.0);
                double fpow = -4*pow(M_PI, 2); 

                double farten;
                double dis;
                double theta;
                double tol = 1e-11;
    
                vector<double> C = {b, h, fpow, hpo, hdiv, 2e30};//a vector with all our constants to save ourselves from many a flop

                p->pos = instal(C, p);
                
    
                ofstream tmpfile;
                tmpfile.open("values3.txt");
       
                for(int ii = 0; ii < n; ii++){
                        p->pos = get_pos(p, C);
                        
                        farten = sqrt(pow(p->pos[3],2)+pow(p->pos[4],2)+pow(p->pos[5],2));
                        dis = sqrt(pow(p->pos[0],2)+pow(p->pos[1],2)+pow(p->pos[2],2));
                        
                        if ((fabs(dis-0.3075000) <= tol) and (fabs(farten-12.44) <= tol)){
	                        theta =  atan(p->pos[1]/p->pos[0]);
		                    tmpfile << 206265*(theta) << " " << h*ii << "\n";
                        }
                }
                tmpfile.close();
                return 0;//if the program runs succesfully, return 0
        }                        
};

