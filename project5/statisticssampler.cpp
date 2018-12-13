#include "system.h"
#include "statisticssampler.h"
#include "lennardjones.h"
#include <iostream>
#include <math.h>

using std::ofstream; using std::cout; using std::endl;

StatisticsSampler::StatisticsSampler()
{

}

void StatisticsSampler::saveToFile(System &system)
{
    // Save the statistical properties for each timestep for plotting etc.
    // First, open the file if it's not open already

    m_file.open("statistics.txt", ofstream::app);
    m_file << system.time() << " " << m_temperature << " " << m_kineticEnergy << " " << m_potentialEnergy << " " << m_kineticEnergy+m_potentialEnergy << " " << m_diffusionConst << " " << m_density <<   " \n";
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
    sampleDiffusionCoef(system);
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

void StatisticsSampler::sampleDiffusionCoef(System &system)
{
    double generalDiffusion = 0;

    for(Atom* atom:system.atoms()){

        double rit2 = sqrt(atom->position.x()*atom->position.x()+atom->position.y()*atom->position.y()+atom->position.z()*atom->position.z());
        double ri02 = sqrt(atom->initialPosition.x()*atom->initialPosition.x()+atom->initialPosition.y()*atom->initialPosition.y()+atom->initialPosition.z()*atom->initialPosition.z());
        //cout << rit2 << " " << ri02 << endl;
        double ri2 = pow(fabs(rit2 - ri02),2);
        generalDiffusion +=  ri2/(6.0*system.time());

    }
    m_diffusionConst = generalDiffusion/system.atoms().size();
}
void StatisticsSampler::sampleDensity(System &system)
{
  double rx=0;
  double ry=0;
  double  rz=0;
  double rxi, ryi, rzi;

  int i =0;
  for (Atom *atomi: system.atoms()){
    int j =0;
    for (Atom *atomj: system.atoms()){
      if (j>i){
	rxi=fabs(atomi->position.x()-atomj->position.x());
	if (rxi>rx)rx=rxi;
	ryi=fabs(atomi->position.y()-atomj->position.y());
	if (ryi>ry)ry=ryi;
	rzi=fabs(atomi->position.z()-atomj->position.z());
	if (rzi>rz)rz=rzi;
      }
      j++;
    }
    i++;
  }
  m_density=system.atoms().size()/(rx*ry*rz);
}
