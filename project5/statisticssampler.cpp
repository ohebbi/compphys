#include "system.h"
#include "statisticssampler.h"
#include "lennardjones.h"
#include <iostream>

using std::ofstream; using std::cout; using std::endl;

StatisticsSampler::StatisticsSampler()
{

}

void StatisticsSampler::saveToFile(System &system)
{
    // Save the statistical properties for each timestep for plotting etc.
    // First, open the file if it's not open already
  
    m_file.open("statistics.txt", ofstream::app);
    m_file << m_temperature << " " << m_kineticEnergy << " " << m_potentialEnergy << " " << m_kineticEnergy+m_potentialEnergy << " \n";
    // Print out values here
    m_file.close();
}

void StatisticsSampler::sample(System &system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    saveToFile(system);
}

void StatisticsSampler::sampleKineticEnergy(System &system)
{
    m_kineticEnergy = 0; // Remember to reset the value from the previous timestep
    for(Atom *atom : system.atoms()) {
      m_kineticEnergy += 0.5*atom->mass()*atom->velocity.lengthSquared();
    }
}

void StatisticsSampler::samplePotentialEnergy(System &system)
{
    m_potentialEnergy = system.potential().potentialEnergy();
}

void StatisticsSampler::sampleTemperature(System &system)
{
    //vec3 systemSize = system.getSystemSize();
    //double N_atoms = 4*systemSize(0)*systemSize(1)*systemSize(2)/latticeConstant;
    m_temperature = (2.0/3.0)*(m_kineticEnergy)/(system.atoms().size());

    // Hint: reuse the kinetic energy that we already calculated
}

void StatisticsSampler::sampleDensity(System &system)
{

}
