#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"

System::System()
{

}

System::~System()
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
  for (Atom *atomi : m_atoms){
    if (atomi->position.x()< -1e-2){ //tolerance for numerical precision
      atomi->position.components[0]+=m_systemSize.x();
      }
    if (atomi->position.x()>= m_systemSize.x()){
      atomi->position.components[0]-=m_systemSize.x();
	}

    if (atomi->position.y()<-1e-2){
      atomi->position.components[1]+=m_systemSize.y();
    }
    if (atomi->position.y()>= m_systemSize.y()){
      atomi->position.components[1]-=m_systemSize.y();
    }

    if (atomi->position.z()<-1e-2){
      atomi->position.components[2]+=m_systemSize.z();
    }
    if (atomi->position.z()>= m_systemSize.z()){
      atomi->position.components[2]-=m_systemSize.z();
    }

    for(Atom *atomj : m_atoms){
      if(atomi!=atomj){
      double dx=atomj->position.x()-atomi->position.x();
      if (dx>m_systemSize.x()*0.5){
	       dx-=m_systemSize.x();
      }
      if (dx<= -m_systemSize.x()*0.5){
	       dx+=m_systemSize.x();
      }
      double dy=atomj->position.y()-atomi->position.y();
      if (dy>m_systemSize.y()*0.5){
	       dy-=m_systemSize.y();
      }
      if (dy<= -m_systemSize.y()*0.5){
	       dy+=m_systemSize.y();
      }
      double dz=atomj->position.z()-atomi->position.z();
      if (dz>m_systemSize.z()*0.5){
	       dz-=m_systemSize.z();
      }
      if (dz<= -m_systemSize.z()*0.5){
	       dz+=m_systemSize.z();
      }
      }
    }
  }
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
    double px=0;
    double py=0;
    double pz=0;
    //vector<double>= momentum;

  for (Atom*atomi:m_atoms){
    px=atomi->mass()*atomi->velocity.x();
    py=atomi->mass()*atomi->velocity.y();
    pz=atomi->mass()*atomi->velocity.z();
    atomi->momentum.set(px,py,pz);

  }
  double p_x=0;
  double p_y=0;
  double p_z=0;
  for (Atom*atomi:m_atoms){
    p_x+=atomi->momentum.x();
    p_y+=atomi->momentum.y();
    p_z+=atomi->momentum.z();
  }
  double x_momentum = p_x/m_atoms.size();
  double y_momentum = p_y/m_atoms.size();
  double z_momentum = p_z/m_atoms.size();
  for (Atom*atomi:m_atoms){
    atomi->momentum.set(atomi->momentum.x()-x_momentum,atomi->momentum.y()-y_momentum,atomi->momentum.z()-z_momentum);
  }
/*
// Check if momentum is preserved
   p_x=0;
   p_y=0;
   p_z=0;
  for (Atom*atomi:m_atoms){
    p_x+=atomi->momentum.x();
    p_y+=atomi->momentum.y();
    p_z+=atomi->momentum.z();
  }
*/

//  std::cout <<  "momentum: "<< m_atoms[i]->momentum << "\n";

}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double b, double temperature) {
    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).
    double b_half = b*0.5;
    double numberOfUnitCells1 = pow(numberOfUnitCellsEachDimension,3);
    //Atom*atom;
    for(int i=0; i<numberOfUnitCellsEachDimension; i++) {
        for(int j=0; j<numberOfUnitCellsEachDimension; j++) {
            for(int k=0; k<numberOfUnitCellsEachDimension; k++) {

                Atom *atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                atom1->position.set(i*b,j*b,k*b);
                atom1->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom1);

                Atom *atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                atom2->position.set(i*b +b_half,j*b+b_half,k*b);
                atom2->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom2);

                Atom *atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                atom3->position.set(i*b,j*b+b_half,k*b+b_half);
                atom3->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom3);

                Atom *atom4 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                atom4->position.set(i*b +b_half,j*b,k*b + b_half);
                atom4->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom4);

}

}

    }
    setSystemSize(vec3(numberOfUnitCellsEachDimension*b, numberOfUnitCellsEachDimension*b, numberOfUnitCellsEachDimension*b)); // Remember to set the correct system size!
}

void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}

void System::step(double dt) {
    m_integrator.integrate(*this, dt);
    m_steps++;
    m_time += dt;
}
