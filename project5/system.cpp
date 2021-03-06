#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"
#include "math.h"
#include <iostream>

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
  for(Atom *atomi:m_atoms){
    for(int i = 0; i <= 3; i++){
      if(atomi->position[i] < -m_systemSize.components[i]*0.5) {
	       atomi->position[i] += m_systemSize.components[i];
      }
      else if(atomi->position[i] >= m_systemSize.components[i]*0.5){
	       atomi->position[i] -= m_systemSize.components[i];
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
    atomi->velocity.set(atomi->momentum.x()/atomi->mass(),atomi->momentum.y()/atomi->mass(),atomi->momentum.z()/atomi->mass());

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


  std::cout <<  "momentum: "<<p_x << " "<< p_y << " " <<p_z << "\n";
*/
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double b, double temperature) {

    double b_half = b*0.5;

    //Atom*atom;
    for(int i=0; i<numberOfUnitCellsEachDimension; i++) {
        for(int j=0; j<numberOfUnitCellsEachDimension; j++) {
            for(int k=0; k<numberOfUnitCellsEachDimension; k++) {

                Atom *atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                atom1->position.set(i*b,j*b,k*b);
                atom1->initialPosition = (atom1->position);
                atom1->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom1);

                Atom *atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                atom2->position.set(i*b +b_half,j*b+b_half,k*b);
                atom2->initialPosition = (atom2->position);
                atom2->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom2);

                Atom *atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                atom3->position.set(i*b,j*b+b_half,k*b+b_half);
                atom3->initialPosition = (atom3->position);
                atom3->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom3);

                Atom *atom4 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                atom4->position.set(i*b +b_half,j*b,k*b + b_half);
                atom4->initialPosition = (atom4->position);
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
