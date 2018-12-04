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
    double systemSize_saveflops = 0.5*systemSize.x(); //ikke generell
    int i=0;
    int j=0;
    for (Atom*atomi:system.atoms()){
      atomi->force.set(0,0,0);
    }
    for (Atom*atomi:system.atoms()){
      double fx =0;
      double fy =0;
      double fz =0;
        for (Atom*atomj:system.atoms()){

            if (atomi != atomj and j>i){
                  double rij = 0;
                  double dr = 0;
                  //Check if the distance between two atoms at different lattices is less
                  //than the distance inside the lattice for the same atoms.
                  double dx=atomj->position.x()-atomi->position.x();
                  double dy=atomj->position.y()-atomi->position.y();
                  double dz=atomj->position.z()-atomi->position.z();

                  if (dx>systemSize_saveflops){
                       dx-=systemSize.x();
                  }
                  if (dx<= -systemSize_saveflops){
                       dx+=systemSize.x();
                  }
                  if (dy>systemSize_saveflops){
                       dy-=systemSize.y();
                  }
                  if (dy<= -systemSize_saveflops){
                       dy+=systemSize.y();
                  }
                  if (dz>systemSize_saveflops){
                       dz-=systemSize.z();
                  }
                  if (dz<= -systemSize_saveflops){
                       dz+=systemSize.z();
                  }


                  dr = sqrt(dx*dx+dy*dy+dz*dz);
                  //std::cout << "dr "<< dr << "\n";
                  rij = sqrt(pow( (atomi-> position.x() - atomj->position.x()), 2) + pow( (atomi-> position.y() - atomj->position.y() ),2) + pow( (atomi-> position.z() - atomj->position.z() ),2));
                  //std::cout << "rij "<< rij << "\n";

                  if (dr < rij) (rij = dr);
                  //std::cout << rij << "da ble det denne" << "\n";

                  double rij8 = pow(rij,8); //kan utbedres
                  double rij14 = pow(rij,14);

                  double dudr = (96/rij14 - 48/rij8);

                  double dfx = dudr*(atomi-> position.x() - atomj->position.x());
                  double dfy = dudr*(atomi-> position.y() - atomj->position.y());
                  double dfz = dudr*(atomi-> position.z() - atomj->position.z());

                  fx+=dfx;
                  fy+=dfy;
                  fz+=dfz;

                  atomj->force.set(atomj->force[0]-dfx,atomj->force[1]-dfy,atomj->force[2]-dfz);

                  //atomj->force.set(-f[0],-f[1],-f[2]);


                  if (std::find(checked.begin(), checked.end(), atomj) == checked.end()){
                    m_potentialEnergy += 4*((rij14/(rij*rij)) - (rij8/(rij*rij)));
                  }

          }
      j += 1;
      checked.push_back(atomi);
      }
    atomi->force.set(fx,fy,fz);
    i +=1;
    }
}
