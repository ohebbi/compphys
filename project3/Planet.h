#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include<fstream>

using namespace std;

//parent class
class Planet {
    public:
        double mass; //mass of the planet
        vector<double> inv; //initial values for position and velocity
        vector<double> pos; //values for position, velocity and acceleration
        Planet(double mass1, vector<double> inv1){
            mass = mass1;
            inv = inv1;
        }
};


//all the planets with inital values gets created here ass subclasses of the parent class Planet

class Sun : public Planet{ //you may think that we believe the sun to be a planet. This is however not the case; for we are on a mission to save time. Therefore, we use the sun as a planet only for practicle purposes. Yours truly; Erlend, Ingvild and Oliver.
    public:            
        Sun():Planet(2e30, {0,0,0,0,0,0}){           
        };   
};

class Mercury : public Planet{
    public:            
        
        Mercury():Planet(3.3e23,{-0.3294, -0.29288, 0.005618, 4.82, -7.13, -1.025}){  //standard         
        };
        /*
        Mercury():Planet(3.3e23,{0.30750000, 0,0,0,12.44000,0}){ //for perhilion          
        };*/
        
            
};

class Venus : public Planet{
    public:            
        Venus():Planet(4.9e24, {0.7244, -0.03279, -0.042427, 0.37, 7.34, 0.079}){           
        };    
};

class Earth : public Planet{
    public:            
        Earth():Planet(6.0e24, {0.980205, 0.205756, -0.0000874989, -1.3497, 6.1340, 0.00042945}){           
        };    
};

class Mars : public Planet{
    public:            
        Mars():Planet(6.6e23, {1.3490045, -0.29756, -0.03956, 1.318, 5.420, 0.081}){   
                   
        };
};

class Jupiter : public Planet{
    public:            
        Jupiter():Planet(1.9e27, {-2.71845, -4.6282558, 0.0800034, 2.342, -1.264, -0.04712}){           
        };
};

class Saturn : public Planet{
    public:            
        Saturn():Planet(5.5e26, {1.497, -9.9435, 0.1133, 1.901, 0.296, -0.0807}){           
        };
};

class Uranus : public Planet{
    public:            
        Uranus():Planet(8.8e25, {17.196, 0.9965, -0.1858, -0.73, 1.175, 0.0139}){           
        };    
};

class Neptune : public Planet{
    public:            
        Neptune():Planet(1.03e26, {28.91, -7.753, -0.5067, 0.289, 1.11, -0.0297}){           
        };    
};

class Pluto : public Planet{
    public:            
        Pluto():Planet(1.31e22, {11.614, -31.58, 0.01979, 1.0975, 0.1535, -0.331}){           
        };   
};










