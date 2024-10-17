#ifdef ISOTROPIC
#error gaussian_patchy_disc.cpp cannot be linked to isotropic VMMC library!
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "System.h"
#include <csignal>
#include <vector>
#include <algorithm>

int main(int argc, char** argv)
{
    std::cout << "# PatchyDisc compiled on " << COMPILED_ON << std::endl;
    std::cout << "# GIT information: branch '"<< GIT_BRANCH <<"', commit hash '" << GIT_COMMIT_HASH << "'" <<std::endl << std::endl;

    if(argc == 1){
        fprintf(stderr, "Usage: %s input.json\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    std::string inputFileName = argv[1];

	//Here we handle a few SIG* signals: whenever we intercept one of these signals
	//the program will quit the main loop and die as gracefully as possible.
	signal(SIGTERM, gbl_terminate);
	signal(SIGABRT, gbl_terminate);
	signal(SIGINT, gbl_terminate);
	signal(SIGUSR2, gbl_terminate);

	//Load the input file and Initialise (inp is a json object) 
    inp = ReadJsonFromFile(inputFileName);
    parseInputFile();
    std::vector<Particle> particles(nParticles);    // particle container
    std::vector<unsigned int> bonds(nParticles ,0); 
    bool isIsotropic[nParticles];                   // whether the potential of each particle is isotropic
    particles.resize(nParticles);                   // Resize particle container.
    Box box(boxSizeVec);
    unsigned int counting = 0;
    std::cout << " the total number of types is " << top.nTypes << std::endl;
    std::cout << " the surface interaction used is " << surface_interaction << " if you need to change, you have to modify gaussian_patchy.cpp and recompile " << std::endl;
    for (unsigned int i=0;i<top.nTypes;i++){
        for (unsigned int j=counting;j<top.Ni[i]+counting;j++){
            particles[j].type=i;
        }
        counting = counting+top.Ni[i];
    }

    printf("# Number of species %7d\n", nTypes);
    printf("# Number of particles %7d\n", nParticles);
    printf("# Number density is %7.4f\n", nParticles/box.Volume);
    printf("%5.4f\n", top.shift[0][0]);
    // Initialise cell list.
    CellList cells;
    cells.setDimension(dimension);
    cells.initialise(box.boxSize, interactionRange);

    // Initialise random number generator.
    MersenneTwister rng;
    rng.setSeed(seed);
    std::cout << "# Random seed set to "<< rng.getSeed()<< std::endl;

    // Initialise input/output class,
    InputOutput io;
    // Create VMD script.
    io.vmdScript(boxSizeVec);
    // Initialise the patchy disc model.
     std::cout << "The interaction is " << interaction << std::endl;
     Model* ppatchyDisc = nullptr;
     if (interaction == "GaussianPatchyDisc")
     {
         GaussianPatchyDisc patchyDisc(box, particles, cells, top, bonds, maxInteractions, interactionEnergy, interactionRange);
         hasFiniteRepulsion = true;
         std::cout << "Gaussian Patchy Disc interaction initialized." << std::endl;
         ppatchyDisc= new GaussianPatchyDisc(box, particles, cells, top, bonds, maxInteractions, interactionEnergy, interactionRange);
     }
     else if (interaction == "KFOpenSurfPatchyDisc")
     {
         KFOpenSurfPatchyDisc patchyDisc(box, particles, cells, top, bonds, maxInteractions, interactionEnergy, interactionRange);
         hasFiniteRepulsion = true;
         std::cout << "Kern Frenkel Patchy Disc interaction with opening and surface interaction initialized." << std::endl;
 //        ppatchyDisc=&patchyDisc;
         ppatchyDisc= new KFOpenSurfPatchyDisc(box, particles, cells, top, bonds, maxInteractions, interactionEnergy, interactionRange); 
     }
     else if (interaction == "GaussianPatchyDiscHR")
     {
         GaussianPatchyDiscHR patchyDisc(box, particles, cells, top, bonds, maxInteractions, interactionEnergy, interactionRange);
         hasFiniteRepulsion = false;
         std::cout << "Gaussian Patchy Disc interaction with hard repulsion initialized." << std::endl;
         ppatchyDisc= new GaussianPatchyDiscHR(box, particles, cells, top, bonds, maxInteractions, interactionEnergy, interactionRange); 
;
     }
     else if (interaction == "GaussianPatchyDiscHRSW")
     {
         GaussianPatchyDiscHRSW patchyDisc(box, particles, cells, top, bonds, maxInteractions, interactionEnergy, interactionRange);
         hasFiniteRepulsion = false;
         std::cout << "Gaussian Patchy Disc HRSW interaction initialized." << std::endl;
         ppatchyDisc= new GaussianPatchyDiscHRSW(box, particles, cells, top, bonds, maxInteractions, interactionEnergy, interactionRange); 

     }
     else
     {
         GaussianPatchyDisc patchyDisc(box, particles, cells, top, bonds, maxInteractions, interactionEnergy, interactionRange);
         ppatchyDisc= new GaussianPatchyDisc(box, particles, cells, top, bonds, maxInteractions, interactionEnergy, interactionRange);
         std::cerr << "[ERROR] Invalid interaction ("<<interaction<<"!\n";
         exit(EXIT_FAILURE);
    }

    // Initialise particle initialisation object
    Initialise initialise;
    if (init_mode == "from_random_conf")
    {
        //Generate a random particle configuration.
        initialise.random(top, particles, cells, box, rng, false);
    }
    else if (init_mode == "from_init_conf")
    {
        initialise.previous(top, particles, cells, box, rng, false);
        loadInitialConfiguration(init_conf, box, particles, cells);
        loadInitialConfigurationPatches(init_patch_conf,particles);
        std::cout<< "# Starting configuration is loaded from '"<< last_conf<<"'. Setting current step to "<< starting_step <<"."<< std::endl;
    }

    // // For debug
    // for (double t=0; t<=360; t+=1)
    // {
    //     double p1[2] = {10., 10.};
    //     double p2[2] = {10.84, 10.84};
    //     double o1[2] = {cos(30./180*M_PI), sin(30./180*M_PI)};
    //     double o2[2] = {cos(t/180*M_PI), sin(t/180*M_PI)};
    //     double E = patchyDisc.computePairEnergy(0, p1, o1, 1, p2, o2);
    //     std::cout << t << " " << E <<"     "<< o2[0]<< " "<<o2[1] << std::endl << std::endl;
    // }
    // exit(1);

    // Initialise data structures needed by the VMMC class.
    double coordinates[dimension*nParticles];
    double orientations[dimension*nParticles];
    // Copy particle coordinates and orientations into C-style arrays.
    for (unsigned int i=0;i<nParticles;i++)
    {
        for (unsigned int j=0;j<dimension;j++)
        {
            coordinates[dimension*i + j] = particles[i].position[j];
            orientations[dimension*i + j] = particles[i].orientation[j];
        }
        // Set all particles as anisotropic.
        isIsotropic[i] = false;
    }

    // Initialise the VMMC callback functions.
    //using namespace std::placeholders;
    //vmmc::CallbackFunctions callbacks;
// #ifndef ISOTROPIC
    //callbacks.energyCallback =
      //  std::bind(&GaussianPatchyDisc::computeEnergy, patchyDisc, _1, _2, _3);
   // callbacks.pairEnergyCallback =
     //   std::bind(&GaussianPatchyDisc::computePairEnergy, patchyDisc, _1, _2, _3, _4, _5, _6);
    // callbacks.interactionsCallback =
     //   std::bind(&GaussianPatchyDisc::computeInteractions, patchyDisc, _1, _2, _3, _4);
    // callbacks.postMoveCallback =
      //  std::bind(&GaussianPatchyDisc::applyPostMoveUpdates, patchyDisc, _1, _2, _3);
// #else
   // callbacks.energyCallback =
    //    std::bind(&GaussianPatchyDisc::computeEnergy, patchyDisc, _1, _2);
    // callbacks.pairEnergyCallback =
     //   std::bind(&GaussianPatchyDisc::computePairEnergy, patchyDisc, _1, _2, _3, _4);
    // callbacks.interactionsCallback =
    //    std::bind(&GaussianPatchyDisc::computeInteractions, patchyDisc, _1, _2, _3);
   // callbacks.postMoveCallback =
     //   std::bind(&GaussianPatchyDisc::applyPostMoveUpdates, patchyDisc, _1, _2);
// #endif

    // Initialise VMMC object.
  //  vmmc::VMMC vmmc(rng, nParticles, dimension, coordinates, orientations,
    //    0.15, 0.2, 0.5, 0.5, maxInteractions, &boxSizeVec[0], isIsotropic, hasFiniteRepulsion, callbacks);
    // Initialise single particle move object.
     SingleParticleMove MC(rng, ppatchyDisc, surface_interaction, 0.4, 0.2, 0.5, 0.75, false);

    // Initialisation depending on the restart_step_counter
    if (restart_step_counter)
    {
        starting_step = 0;
        bool clearFile = true;
        std::cout<< "# Resetting the step counter to zero. Owerwiting '"<< trajectory << "' and '"<< log_file<<"'.\n" << std::endl;
        io.appendXyzTrajectory(trajectory, patchfile, starting_step, box, particles, top, clearFile, true);
        io.appendLog(log_file, starting_step, ppatchyDisc->getEnergy(), clearFile);
    }
    else
    {
        bool clearFile = false;
        std::cout<< "# Resuming the simulation; Appending to '"<< trajectory << "' and '"<< log_file<<"'.\n"<< std::endl;
        if (starting_step % output_every == 0)
        {
            io.appendXyzTrajectory(trajectory, patchfile, starting_step, box, particles, top, clearFile, true);
            io.appendLog(log_file, starting_step, ppatchyDisc->getEnergy(), clearFile);
        }
    }

    // Execute the simulation.
    for(curr_step = starting_step; curr_step < sweeps && !stop; curr_step++)
    {
        for (unsigned int i=0; i<nParticles; i++)
        {
           // vmmc ++;
            MC ++;   // both moves cannot be used at the same time. BUG?
        }
        if(curr_step > 0 && (curr_step % (output_every)) == 0)
        {
            io.appendXyzTrajectory(trajectory, patchfile, curr_step, box, particles, top, false, true);
            io.appendLog(log_file, curr_step, ppatchyDisc->getEnergy(), false);
            io.appendPatchstate(last_patch_states_file,particles,top,true);
        }
    }
    printf("\nWriting the last conf to `%s`.", last_conf.c_str());
    io.appendXyzTrajectory(last_conf, last_patchfile, curr_step, box, particles, top, true, true);
    std::cout << "\nComplete!\n";

    // We're done!
    return (EXIT_SUCCESS);
}
