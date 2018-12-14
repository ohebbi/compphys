#include "lennardjones.h"
#include "system.h"
#include<math.h>
#include<algorithm>
#include<iostream>

using namespace std;

double LennardJones::potentialEnergy() const
{
    return m_potentialEnergy;
}

double LennardJones::sigma() const
{
    return m_sigma;
}

void LennardJones::setSigma(double sigma)
{
    m_sigma = sigma;
}

double LennardJones::epsilon() const
{
    return m_epsilon;
}

void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void LennardJones::calculateForces(System &system)
{
    m_potentialEnergy = 0;

    int i = 0;

    for(Atom *atomi:system.atoms()) {
      int j = 0;

      for(Atom *atomj:system.atoms()) {
 	          if (j > i){

                vec3 rij = atomi->position - atomj->position;

	              for (int k = 0; k < 3; k++){
	                   if (rij[k] <= - system.systemSize()[k]*0.5){
                          rij[k] += system.systemSize()[k];
                      }
                      else if (rij[k] >= system.systemSize()[k]*0.5){
                          rij[k] -= system.systemSize()[k];
                      }
                }
            double rij2 = rij.lengthSquared();
            if (rij2 <= 9){
                double rij2inv = 1.0/rij2; // inverse delta r, unitless, squared
                double rij2inv3 = pow(rij2inv, 3); //attractive
                double rij2inv6 = pow(rij2inv3, 2); //repulsive

                m_potentialEnergy += 4*(rij2inv6 - rij2inv3);

                vec3 F = rij/rij2*(48*rij2inv6 - 24*rij2inv3);

                atomj->force -= F;
                atomi->force += F;
            }
	         }
        j++;

        }
      i++;
    }
}
