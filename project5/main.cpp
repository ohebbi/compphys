#include "math/random.h"
#include "lennardjones.h"
#include "velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>

using namespace std;

int main(int numberOfArguments, char **argumentList)
{

    for(int temp = 500; temp <= 3000; temp+=25)
    {
    cout << temp << endl;
    ofstream m_file;
    m_file.open("statistics.txt");
    m_file.close();
    int numberOfUnitCells =5;
    double initialTemperature = UnitConverter::temperatureFromSI(temp); // measured in Kelvin
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms


    // If a first argument is provided, it is the number of unit cells
    if(numberOfArguments > 1) numberOfUnitCells = atoi(argumentList[1]);
    // If a second argument is provided, it is the initial temperature (measured in kelvin)
    if(numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));
    // If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
    if(numberOfArguments > 3) latticeConstant = UnitConverter::lengthFromAngstroms(atof(argumentList[3]));

    double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds.

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    cout << "One unit of diffusion coeffisient is " << UnitConverter::diffusionToSI(1.0) << "" << endl;

    System system;
    system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature);
    system.potential().setEpsilon(1.0/8.6173303e-5); // /8.6173303e-5
    system.potential().setSigma(1.0);


    system.removeTotalMomentum();

    StatisticsSampler statisticsSampler;
    IO movie("movie.xyz"); // To write the state to file

    cout << setw(20) << "Timestep" <<
            setw(20) << "Time" <<
            setw(20) << "Temperature" <<
            setw(20) << "KineticEnergy" <<
            setw(20) << "PotentialEnergy" <<
            setw(20) << "TotalEnergy" <<
            setw(20) << "Diffusionconstant" << endl;

    for(int timestep=0; timestep<10000; timestep++) {
        statisticsSampler.sample(system);
        system.step(dt);
        if(system.steps() % 1000 == 0 ) {
          // Print the timestep every 100 timesteps
          cout << setw(20) << system.steps() <<
            setw(20) << system.time() <<
            setw(20) << statisticsSampler.temperature() <<
            setw(20) << statisticsSampler.kineticEnergy() <<
            setw(20) << statisticsSampler.potentialEnergy() <<
            setw(20) << statisticsSampler.totalEnergy() <<
            setw(20) << statisticsSampler.diffusionConst() << endl;

        }
        movie.saveState(system);
    }

    movie.close();
    ifstream stream1;
    stream1.open("statistics.txt", ios::in);
    string line;
    ofstream stream2;
    stream2.open(to_string(temp) + "diffusion.txt", ios::out);
    while(getline(stream1,line)){
        stream2 << line << " \n";
    }
    stream1.close();
    stream2.close();

    }
    return 0;
}
