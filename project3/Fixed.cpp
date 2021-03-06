#include "Planet.h"
#include "Solver_fixed.h"
#include "Solver_fixed_perhilion.h"

int main(int argc,char* argv[]){

    double final_time = 30; //in years

    double betta = 3.0; //for varying the gravitational force in task D
    int n = 1e6; //number of integration points

    //creating the planets from Planet.h
    Sun* sol = new Sun();
    Mercury* mercury = new Mercury();
    Venus* venus = new Venus();
    Earth* earth = new Earth();
    Mars* mars = new Mars();
    Jupiter* jupiter = new Jupiter();
    Saturn* saturn = new Saturn();
    Uranus* uranus = new Uranus();
    Neptune* neptune = new Neptune();
    Pluto* pluto = new Pluto();   
    
    vector<Planet*> planets = {mercury, venus, earth, mars, jupiter, saturn, uranus, neptune, pluto}; //vector containing all the objects you want to use in the simulation
           
    Solver* solver = new Solver(); //creates a solver
    solver->solution(final_time, betta, planets, n);//solves the differential equations for the given system

}
