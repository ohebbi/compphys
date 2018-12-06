#include "lennardjones.h"
#include "system.h"
#include<math.h>
#include<algorithm>
#include<iostream>


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
    
    system.applyPeriodicBoundaryConditions();
    for(Atom *atom : system.atoms()) {
        atom->resetForce();
    }
    
    std::vector<Atom*> checked;
    
    vec3 systemSize = system.getSystemSize();
    
   
    int i = 0;
    double dfx;
    double dfy;
    double dfz;

    m_potentialEnergy=0;
    for (Atom*atomi:system.atoms()){

      double fx =0;
      double fy =0;
      double fz =0;
      int j = 0;

        for (Atom*atomj:system.atoms()){
            
                
            if (atomi->position.x() != atomj->position.x() and atomi->position.y() != atomj->position.y() and atomi->position.z() != atomj->position.z()){
                double rij = 0;

                //Check if the distance between two atoms at different lattices is less
                //than the distance inside the lattice for the same atoms.
                  
                double dx = (atomj-> position.x() - atomi->position.x());
                if(dx > systemSize.x()*0.5){
                    dx -= systemSize.x();
                }
                 if(dx <= -systemSize.x()*0.5){
                    dx += systemSize.x();
                }
                double dy = (atomj-> position.y() - atomi->position.y());
                if(dy > systemSize.y()*0.5){
                    dy -= systemSize.y();
                }
                 if(dy <= -systemSize.y()*0.5){
                    dy += systemSize.y();
                }
                double dz = (atomj-> position.z() - atomi->position.z());
                if(dz > systemSize.z()*0.5){
                    dz -= systemSize.z();
                }
                 if(dz <= -systemSize.z()*0.5){
                    dz += systemSize.z();
                }                       
                                
                rij = sqrt(dx*dx+dy*dy+dz*dz);
                
                double rij8 = pow(rij,8); //kan utbedres
                double rij14 = pow(rij,14); 
                double sigma12 = pow(m_sigma, 12);
                double sigma6 = pow(m_sigma, 6);
                
                          
                      
                double dudr = (48.0*sigma12/rij14 - 24.0*sigma6/rij8);
                if(fabs(dx)< 2.5 and fabs(dx) > 1.0){
                  dfx = dudr*dx;
                    }
                else{
                    //std::cout << dudr*dx << std::endl;
                    dfx = 0;

                }
                if(fabs(dy)<2.5 and fabs(dy) > 1.0){
                   dfy = dudr*dy;
                    }
                else{
                    //std::cout << dudr*dy << std::endl;
                    dfy = 0;
                    
                }
                if(fabs(dz)<2.5 and fabs(dz) > 1.0){
                 dfz =dudr*dz;                          
                    }
                else{
                    //std::cout << dudr*dz << std::endl;
                    dfz = 0;
                }

                if(fabs(dfx) > 100 or fabs(dfy) > 100 or fabs(dfz) > 100){
                    std::cout << dudr << std::endl;
                }
                fx+=dfx;
                fy+=dfy;
                fz+=dfz;
                //atomj->force.set(atomj->force.components[0]-dfx,atomj->force.components[1]-dfy,atomj->force.components[2]-dfz);
                
                 
                  


                  if (std::find(checked.begin(), checked.end(), atomj) == checked.end()){
                    m_potentialEnergy += 4*m_epsilon*(((rij*rij)/rij14) - ((rij*rij)/rij8));
                  }

        
            } 
     
      
        j += 1;
        }
        i += 1; 
    checked.push_back(atomi);
    vec3 allForces = {fx,fy,fz};

    atomi->force = m_epsilon*allForces;
    }
}
