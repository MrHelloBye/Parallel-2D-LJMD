//This is a parallelized (using MPI) object oriented molecular dynamics code which relies on many classes. This file allows to control most of the
//program parameters and choose which simulation approach (NVT, NVE) to use. Certain blocks need
//to be commented/uncommented to use different approaches.

//Output is a movie.xyz file with atom positions which can allow to create animations in i.e Ovito.
//Additional system properties such as energies, temperature, pressure, density, diffusion coefficient can be output if desired
//by using StatisticsSampler class. Output can be made in SI of Lennard Jones units by using UnitConverter class.
//The frequency of output can be varied as desired. Of course, less frequent output, makes the program run faster.

//Author: Timofey Golubev with help of Xukun Xiang based on the code found at https://github.com/andeplane/molecular-dynamics-fys3150.

#include "global.h"
#include "random.h"
#include "lennardjones.h"
#include "velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <time.h>
#include <stdio.h> //for printf
#include <mpi.h>  //get the mpi functions from here

#include "../comms.hpp"
#include "../Controller/controllerState.hpp"


using namespace std;
using namespace chrono;

int modelMain(int argc, char **argv)
{
    int my_id, nprocs;
    
    //each processor finds out what it's ID (aka rank) is and how many processors there are
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);   //find ID
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);  //find # of processors

    double initialTemperature = 600.;//in K
    double currentTemperature;
    double latticeConstant =3.8;    //in angstroms
    double sigma = 3.4;             //atom/particle diameter in Angstrom for LJ potential
    double epsilon = 1.0318e-2;     // epsilon from LJ in eV
    double side_length;

    double total_dt_time= 0.0;

    //all variables will be defined in EACH processor
    vec2 Total_systemSize(90,30); //since using SC lattice--> just gives # of atoms in each dimension
    vec2 subsystemSize;
    subsystemSize[0] = Total_systemSize[0]/(nprocs-1);  //1D parallelization along x
    subsystemSize[1] = Total_systemSize[1]; //MUST CHANGE THIS TO /nprocs when do 2D parallelization

    vec2 subsystemOrigin; //bottom left corner of each subsystem, this will be defined in each processor seperately

    int StatSample_freq = 10;

    int N_time_steps = 10000; //number of time steps

    //for NVT ensemble
    int N_steps = 1;  //number of steps over which to gradually rescale velocities: has to be large enough to prevent instability
    //using the more basic thermostat of sharp rescaling every 100 steps

    //-----------------------------------------------------------------------------------------------------------------------------------------------
    //Initialize MD units
    UnitConverter::initializeMDUnits(sigma, epsilon);
    //double mass = UnitConverter::massFromSI(6.63352088e-26); // mass in kg
    initialTemperature = UnitConverter::temperatureFromSI(initialTemperature);
    latticeConstant = UnitConverter::lengthFromAngstroms(latticeConstant);
    double dt = UnitConverter::timeFromSI(2e-15); //

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    cout <<"Epsilon is " << epsilon <<" eV and sigma is " <<sigma <<" in Angstrom" << endl;
    cout << "mass in LJ units" << mass <<endl;
    cout << "timestep in LJ units" <<dt <<endl;

    System system;
    if(my_id !=0){
        system.createSCLattice(Total_systemSize,subsystemSize, latticeConstant, initialTemperature, mass, subsystemOrigin);
        system.potential().setEpsilon(1.0); //if don't set to 1.0, must modify LJ eqn.
        system.potential().setSigma(UnitConverter::lengthFromAngstroms(sigma));      //i.e. LJ atom/particle diameter,
        system.m_sample_freq=100; //statistics sampler freq.
        system.removeTotalMomentum();
    }

    //create_MPI_ATOM();

    StatisticsSampler statisticsSampler;

    FILE *movie;
    //IO movie("movie.xyz");

    if(my_id != 0){
        cout << setw(20) << "Timestep" <<
                setw(20) << "Time" <<
                setw(20) << "Temperature(not K!)" <<
                setw(20) << "KineticEnergy" <<

                setw(20) << "PotentialEnergy" <<
                setw(20) << "TotalEnergy" << endl;
    }else{
        //movie = fopen("movie.xyz","w+");
        //NumInBox_array = new int[nprocs]; //allocate memery only in proc 1, array for recieving number of particles in each procs domain
    }

    high_resolution_clock::time_point start2 = high_resolution_clock::now();

    ControllerState controllerState;

    for(int timestep=0; timestep< N_time_steps; timestep++) {
        if(my_id != 0){
            system.step(dt);  //only do timestepping for non-root procs


            //use sampler to calculate system parameters periodically
            if(timestep % system.m_sample_freq ==0){
                statisticsSampler.sample(system);
            }


            //periodically rescale Velocities to keep T constant (NVT ensemble)
            if(timestep % 100 == 0){
                statisticsSampler.sampleTemperature(system); //sample temperature, so can rescale as often as we want
                currentTemperature = statisticsSampler.temperature();  //this gets the value of temperature member variable
                //Note: initial temperature is the desired temperature here
                system.rescaleVelocities(statisticsSampler, currentTemperature, initialTemperature, N_steps);
            }

            if( timestep % 100 == 0 ) {
                // Print the timestep and system properties every 100 timesteps
                /*cout << setw(20) << system.steps() <<
                        setw(20) << system.time() <<
                        setw(20) << statisticsSampler.temperature() <<
                        setw(20) << statisticsSampler.kineticEnergy() <<
                        setw(20) << "proc" << my_id <<
                        setw(20) << statisticsSampler.potentialEnergy() <<
                        setw(20) << statisticsSampler.totalEnergy() << endl;*/
                cout<<"Model: " <<my_id<<" Controller: "<<controllerState<<endl;
            }

            //Send the model data to the view (root)
            comms_sendModelData(system);

            //Get the controller state
            comms_recvControllerState(&controllerState);
        }

    }//end of for loop

    //timing
    high_resolution_clock::time_point finish2 = high_resolution_clock::now();
    duration<double> time2 = duration_cast<duration<double>>(finish2-start2);
    cout << "Total CPU time for proc " << my_id <<" is " << time2.count() << endl;

    MPI_Finalize();

    return 0;
}
