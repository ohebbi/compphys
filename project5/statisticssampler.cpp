#include "system.h"
#include "statisticssampler.h"
#include "lennardjones.h"
#include <iostream>
#include <math.h>

using std::ofstream; using std::cout; using std::endl;

StatisticsSampler::StatisticsSampler()
{

}

void StatisticsSampler::saveToFile(System &system, double T)
{
    // Save the statistical properties for each timestep for plotting etc.
    // First, open the file if it's not open already
    m_file.open("statistics-diff.txt", ofstream::app);
    m_file << system.time() << " " << m_temperature << " " << m_kineticEnergy << " " << m_potentialEnergy << " " << m_kineticEnergy+m_potentialEnergy << " " << m_diffusionConstant << "\n";
    // Print out values here
    m_file.close();
}

void StatisticsSampler::sample(System &system, double T)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    sampleDiffusionConstant(system);
    saveToFile(system, T);
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
    m_temperature = (2.0/3.0)*(m_kineticEnergy)/(system.atoms().size());

}

void StatisticsSampler::sampleDensity(System &system)
{

}

void StatisticsSampler::sampleDiffusionConstant(System &system)
{
  double rit2 = 0;
  double ri02 =0;
  for (Atom *atom: system.atoms()){
      rit2 += atom->position.x()*atom->position.x()+atom->position.y()*atom->position.y()+atom->position.z()*atom->position.z();
      ri02 += atom->initialPosition.x()*atom->initialPosition.x()-atom->initialPosition.y()*atom->initialPosition.y()+atom->initialPosition.z()*atom->initialPosition.z();
  }
  double ri2 = (rit2 - 2*sqrt(rit2*ri02) + ri02)/system.atoms().size();
  m_diffusionConstant = ri2/(6*system.time()*1.00224e3); //*1.00224e3 to scale to cm^2/s
  //cout << m_diffusionConstant << endl;
}
