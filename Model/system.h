#ifndef SYSTEM_H
#define SYSTEM_H
#include "atom.h"
#include "vec2.h"
#include <vector>
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"

class System
{
private:
    //m_ stands for member
    vec2 m_systemSize;
    //vec2 m_subsystemSize;
    vec2 m_halfsystemSize;
    VelocityVerlet m_integrator;
    std::vector<Atom*> m_atoms;  //vector of atom pointers
    LennardJones m_potential;
    double m_time = 0;
    int m_steps = 0;
    //int m_sample_freq;  //no need to make this private b/c actually want to be able to change it

public:
    System();
    ~System();
    int m_sample_freq;
    //void createFCCLattice(vec2 Total_systemSize, double latticeConstant, double temperature,  double mass);
    void createSCLattice(vec2 Total_systemSize, vec2 subsystemSize, double latticeConstant, double temperature,  double mass, vec2 subsystemOrigin);
    void createRandomPositions(int num_particles, double side_length, double temperature,  double mass);
    void applyPeriodicBoundaryConditions();
    void rescaleVelocities(StatisticsSampler &statisticsSampler, double currentTemperature, double desiredTemperature, int N_steps);
    void removeTotalMomentum();
    void removeEscapedAtoms();
    void increaseTemperature(StatisticsSampler &statisticsSampler, double increment);
    void calculateForces();
    void step(double dt);

    // Setters and getters
    std::vector<Atom *> &atoms() { return m_atoms; } // Calling atoms(): Returns a reference to the std::vector of atom pointers
    Atom* &atoms(int index) {return m_atoms[index];}
    double volume() { return m_systemSize[0]*m_systemSize[1]*m_systemSize[2]; }
    vec2 systemSize() { return m_systemSize; }
    //vec2 subsystemSize() { return m_subsystemSize; }
    double systemSize(int j) { return m_systemSize[j]; }
    //double subsystemSize(int j) { return m_subsystemSize[j]; }
    double halfsystemSize(int j){return m_halfsystemSize[j];}
    void setSystemSize(vec2 systemSize) {
            m_systemSize = systemSize;
            m_halfsystemSize = 0.5*systemSize;
        }
    /*
   void setSubSystemSize(vec2 subsystemSize) {
          m_subsystemSize = subsystemSize;
        }
        */
    LennardJones &potential() { return m_potential; }
    double time() { return m_time; }
    void setTime(double time) { m_time = time; }
    VelocityVerlet &integrator() { return m_integrator; }
    int steps() { return m_steps; }
    void setSteps(int steps) { m_steps = steps; }
    int num_atoms() {return m_atoms.size();}
    //void setSampleFreq(int frequency){m_sample_freq = frequency;}
    int sample_freq(){return m_sample_freq;}
};
#endif
