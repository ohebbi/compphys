#include "lennardjones.h"
#include "system.h"
#include<math.h>
#include<algorithm>

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
    std::vector<Atom*> checked;
    vec3 systemSize = system.getSystemSize();
    for (Atom*atomi:system.atoms()){
        for (Atom*atomj:system.atoms()){
            if (atomi != atomj){
                  //Check if the distance between two atoms at different lattices is less
                  //than the distance inside the lattice for the same atoms.
                  double dx=atomj->position.x()-atomi->position.x();
                  double dy=atomj->position.y()-atomi->position.y();
                  double dz=atomj->position.z()-atomi->position.z();

                  if (dx>systemSize.x()*0.5){
                       dx-=systemSize.x();
                  }
                  if (dx<= -systemSize.x()*0.5){
                       dx+=systemSize.x();
                  }
                  if (dy>systemSize.y()*0.5){
                       dy-=systemSize.y();
                  }
                  if (dy<= -systemSize.y()*0.5){
                       dy+=systemSize.y();
                  }
                  if (dz>systemSize.z()*0.5){
                       dz-=systemSize.z();
                  }
                  if (dz<= -systemSize.z()*0.5){
                       dz+=systemSize.z();
                  }


                  double dr = sqrt(dx*dx+dy*dy+dz*dz);
                  std::cout << "dr "<< dr << "\n";
                  double rij = sqrt(pow( (atomi-> position.x() - atomj->position.x()), 2) + pow( (atomi-> position.y() - atomj->position.y() ),2) + pow( (atomi-> position.z() - atomj->position.z() ),2));
                  std::cout << "rij "<< rij << "\n";

                  if (dr < rij) (rij = dr);
                  std::cout << rij << "da ble det denne" << "\n";

                  double rij8 = pow(rij,8); //kan utbedres
                  double rij14 = pow(rij,14);

                  double dudr = (96/rij14 - 48/rij8);

                  double fx=dudr*(atomi-> position.x() - atomj->position.x());
                  double fy=dudr*(atomi-> position.y() - atomj->position.y());
                  double fz=dudr*(atomi-> position.z() - atomj->position.z());
                  atomi->force.set(fx,fy,fz);

                  if (std::find(checked.begin(), checked.end(), atomj) == checked.end()){
                    m_potentialEnergy += 4*((rij14/(rij*rij)) - (rij8/(rij*rij)));
                  }

                   // Remember to compute this in the loop //done

          }
      checked.push_back(atomi);
      }

    }
}
