#include "Planet.h"
#include "Solver.h"

int main(int argc,char* argv[]){
    double final_time = 250;
    double b = 3.0;
    int n = 1e6;
    
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
    
    vector<Planet*> planets = {sol, mercury, venus, earth, mars, jupiter, saturn, uranus, neptune, pluto};
           
    Solver* solver = new Solver();
    solver->bane(final_time, b, planets, n);

}