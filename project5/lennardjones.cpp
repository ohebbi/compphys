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
	  double dx =0;
	  double dy =0;
	  double dz =0;            
                
	  if (atomi!= atomj){
                double rij = 0;

                //Check if the distance between two atoms at different lattices is less
                //than the distance inside the lattice for the same atoms.
                  
                dx = (atomj-> position.x() - atomi->position.x());
		dy = (atomj-> position.y() - atomi->position.y());
		dz = (atomj-> position.z() - atomi->position.z());
		/*
		dx = fabs(dx);
		dx -= static_cast<int>(dx*(1.0/systemSize.x()+0.5)*systemSize.x());
		*/
                
		
		if(dx > systemSize.x()*0.5){
                    dx -= systemSize.x();
                }
                 if(dx <= -systemSize.x()*0.5){
                    dx += systemSize.x();
                }
          
                if(dy > systemSize.y()*0.5){
                    dy -= systemSize.y();
                }
                 if(dy <= -systemSize.y()*0.5){
                    dy += systemSize.y();
                }
                
                if(dz > systemSize.z()*0.5){
                    dz -= systemSize.z();
                }
                 if(dz <= -systemSize.z()*0.5){
                    dz += systemSize.z();
                }
		
                                
                rij = sqrt(dx*dx+dy*dy+dz*dz);
                
                double rij8 = pow(rij,8); //kan utbedres
                double rij14 = pow(rij,14); 
        
                
                
                          
                      
                double dudr = (48.0/rij14-24.0/rij8);
                if(rij < 5.0){
                  dfx = dudr*dx;
                    }
                else{
                    //std::cout << dudr*dx << std::endl;
                    dfx = 0;

                }
                if(rij < 5.0){
                   dfy = dudr*dy;
                    }
                else{
                    //std::cout << dudr*dy << std::endl;
                    dfy = 0;
                    
                }
                if(rij < 5.0){
                 dfz =dudr*dz;                          
                    }
                else{
                    //std::cout << dudr*dz << std::endl;
                    dfz = 0;
                }

                
                fx+=dfx;
                fy+=dfy;
                fz+=dfz;
                //atomj->force.set(atomj->force.components[0]-dfx,atomj->force.components[1]-dfy,atomj->force.components[2]-dfz);
                           
                  if (std::find(checked.begin(), checked.end(), atomj) == checked.end()){
                    m_potentialEnergy += 4*(((rij*rij)/rij14) - ((rij*rij)/rij8));
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
