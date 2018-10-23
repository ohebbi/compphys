#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include<fstream>

using namespace std;

//this function finds the center off mass of the system, and then returns the initial values for the sun
vector<double> find_center(vector<Planet*> planets){
        
        vector<double> center {0,0,0,0,0,0};
        double tot_mass = 0;

        for (int i = 0; i < planets.size(); i++){

                for (int j = 0; j < 3; j++){                
                    center[j] -= planets[i]->mass*planets[i]->inv[j];
                }
                for (int jj = 3; jj < center.size(); jj++){
                    center[jj] -= planets[i]->mass*planets[i]->inv[jj]/planets[0]->mass;
                }
                tot_mass += planets[i]->mass;
        }        
        for (int ii = 0; ii < 3;ii++){
                center[ii] /= tot_mass;
        }
        return center;
}

//here the decomposed force on a planet is calculated
vector<double> force(vector<Planet*> planets, double msun){
        double rej;
        vector<double> f = {0,0,0};
        for(int i = 1; i < planets.size(); i++){
                rej = pow(sqrt(pow((planets[i]->pos[0]-planets[0]->pos[0]),2)+pow((planets[i]->pos[1]-planets[0]->pos[1]),2)+pow((planets[i]->pos[2]-planets[0]->pos[2]),2)),3);

                for(int j = 0; j<f.size(); j++){
                        f[j] += (planets[0]->pos[j]-planets[i]->pos[j])*planets[i]->mass/(rej*msun);
                }
        }
        return f;
}

//instal creates the initial acceleration needed for the Verlet method
vector<double> instal(vector<double> C, vector<Planet*> planets, int in){

        swap(planets[0], planets[in]);//this is important to include all the planets
        
        Planet* first = planets[0];

        for (int i = 0; i< planets.size(); i++){
                planets[i]->pos = planets[i]->inv;
        }
        vector<double> f = force(planets, C[5]);

        double ax1 = C[2]*(f[0]);
        double ay1 = C[2]*(f[1]);
        double az1 = C[2]*(f[2]);

        return {first->pos[0], first->pos[1], first->pos[2], first->pos[3], first->pos[4], first->pos[5], ax1, ay1, az1};
}

//this function returns both position, velocity and acceleration for the given planet
vector<double> get_pos(vector<Planet*> planets, int in, vector<double> C){
        swap(planets[0], planets[in]);//this is important to include all the planets

        Planet* p = planets[0];
        
        for(int i = 0; i < 3; i++){
                p->pos[i] += C[1]*p->pos[i+3]+C[3]*p->pos[i+6];
        }
        
        vector<double> f = force(planets, C[5]);

	double ax = C[2]*f[0];
        p->pos[3] += C[4]*(ax+p->pos[6]);

        double ay = C[2]*f[1];
        p->pos[4] += C[4]*(ay+p->pos[7]);

        double az = C[2]*f[2];
        p->pos[5] += C[4]*(az+p->pos[8]);

        return {p->pos[0], p->pos[1], p->pos[2], p->pos[3], p->pos[4], p->pos[5], ax, ay, az};
}

//this the class that will solve our differential equation
class Solver{
    
    public:
            int solution(float final_time, double b, vector<Planet*> planets, int n){
                planets[0]->inv = find_center(planets);//giving the sun initial values    
            
                double h = final_time/n;     
                double hpo = pow(h,2.0)/2.0;
                double hdiv = (h/2.0);
                double fpow = -4*pow(M_PI, 2); 
    
                vector<double> C = {b, h, fpow, hpo, hdiv, 2e30};//a vector with all our constants to save ourselves from many a flop

                for(int i = 0; i < planets.size(); i++){
                        planets[i]->pos = instal(C, planets, i);
                }
    
                ofstream tmpfile;
                tmpfile.open("values3.txt");
       
                for(int ii = 0; ii < n; ii++){
                        for (int j = 0; j < planets.size(); j++){
                                planets[j]->pos = get_pos(planets, j, C);
                        }
                        if(ii%1000==0){  
                                for (int jj = 0; jj < planets.size(); jj++){
                                        tmpfile <<  planets[jj]->pos[0] << " " << planets[jj]->pos[1] << " " << planets[jj]->pos[2] << " " << jj << " " << "\n";//printing the position of all the planets to a text file that shall be read by a python program
                                }
 	                }
                }
                tmpfile.close();
                return 0;//if the program runs succesfully, return 0
        }                        
};

