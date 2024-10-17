#ifndef _SYSTEM_H
#define _SYSTEM_H

#include "Box.h"
#include "CellList.h"
#include "Initialise.h"
#include "InputOutput.h"
#include "Model.h"
#include "Particle.h"
#include "Top.h"
#include <csignal>
#include "Utils.h"

#include "SingleParticleMove.h"
#include "VMMC.h"

#include "PotentialsList.h"
#include <random>
#include <vector>

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif

#ifndef GIT_COMMIT_HASH
#define GIT_COMMIT_HASH "?"
#endif
#ifndef GIT_BRANCH
#define GIT_BRANCH "?"
#endif
#ifndef COMPILED_ON
#define COMPILED_ON "?"
#endif


// FORWARD DECLARATIONS

class  Box;
class  CellList;
struct Particle;


int stop = 0;
void gbl_terminate(int arg) {
	/**
	 * The next time the signal is intercepted, the default handler will be invoked
	 * in order to avoid making the program unkillable.
	 */
	signal(arg, SIG_DFL);
	fprintf(stderr, "# Caught SIGNAL %d; setting stop = 1\n", arg);
	stop = 1;
}

// Simulation parameters.
unsigned int dimension = 2;                     // dimension of simulation box
double interactionRange = 1.0;                  // interaction range. #diameter of patch (in units of particle diameter)
double interactionEnergy = 1.0;                 // pair interaction energy scale (in units of kBT)
unsigned int maxInteractions = 36;              // maximum number of interactions per particle


nlohmann::json inp;                             // inputfile json
unsigned int nParticles;
unsigned int nTypes;
double surface_interaction=0; // Right now surface interaction is not 
Top top;
std::vector<double> boxSizeVec (dimension);
std::string init_conf="init_conf.xyz";
std::string init_patch_conf="init_patch_conf.xvg";
std::string trajectory = "trajectory.xyz";
std::string patchfile = "patch_status.xvg";
std::string last_patch_states_file="last_states.xvg";
std::string last_patchfile = "last_patch.xvg";
std::string last_conf = "last_conf.xyz";
std::string init_mode = "from_init_conf";
std::string log_file = "log.dat";

unsigned int output_every = 1000;
unsigned int sweeps = 10000;
long long int starting_step = 0;
long long int curr_step = 0;
bool restart_step_counter=true;

std::string interaction = "GaussianPatchyDisc";
bool hasFiniteRepulsion;

unsigned int seed = std::random_device{}();

void getParamsFromInitConf(std::string fileName)
{
    std::ifstream dataFile;
    std::vector<double> dummy (dimension);
    unsigned int type, max_type=0;
    dataFile.open(fileName.c_str());

    // Check that the file is valid.
    if (dataFile.good())
    {
        dataFile >> nParticles;
        dataFile >> starting_step;
        for (unsigned int j=0;j<dimension;j++) dataFile >> boxSizeVec[j];
        // find out nTypes
        for (unsigned int i=0;i<nParticles;i++)
        {
            dataFile >> type;
            if (type>max_type) max_type=type;
            for (unsigned int j=0;j<dimension;j++) dataFile >> dummy[j];
            for (unsigned int j=0;j<dimension;j++) dataFile >> dummy[j];
        }
        nTypes = max_type + 1;
    }
    else
    {
        std::cerr << "[ERROR] getParamsFromInitConf(): Invalid init_conf file ('"<< fileName <<"')!\n";
        exit(EXIT_FAILURE);
    }
    // Close file stream.
    dataFile.close();
}

void getParamsFromInitConf2(std::string fileName)
{
    std::ifstream dataFile;
    std::vector<double> dummy (dimension);
    unsigned int type;
    dataFile.open(fileName.c_str());

    // Check that the file is valid.
    if (dataFile.good())
    {
        dataFile >> nParticles;
        dataFile >> starting_step;
        for (unsigned int j=0;j<dimension;j++) dataFile >> boxSizeVec[j];
        // find out nTypes
        for (unsigned int i=0;i<nParticles;i++)
        {
            dataFile >> type;
            top.Ni[type]+=1;
            for (unsigned int j=0;j<dimension;j++) dataFile >> dummy[j];
            for (unsigned int j=0;j<dimension;j++) dataFile >> dummy[j];
        }
    }
    else
    {
        std::cerr << "[ERROR] getParamsFromInitConf(): Invalid init_conf file ('"<< fileName <<"')!\n";
        exit(EXIT_FAILURE);
    }
    // Close file stream.
    dataFile.close();
}


void loadInitialConfiguration(std::string fileName, Box& box, std::vector<Particle>& particles, CellList& cells)
{
    std::ifstream dataFile;
    std::vector<double> dummy (dimension);
    for (unsigned int i=0; i<nTypes; i++) top.Ni[i]=0;

    dataFile.open(fileName.c_str());

    // Check that the file is valid.
    if (dataFile.good())
    {
        dataFile >> nParticles;
        dataFile >> starting_step;
        for (unsigned int j=0;j<dimension;j++) dataFile >> boxSizeVec[j];
        for (unsigned int i=0;i<nParticles;i++)
        {
            // Set particle index.
            particles[i].index = i;

            // Set particle index, and increment the particle number container.
            dataFile >> particles[i].type;
            top.Ni[particles[i].type] += 1;

            // Resize position and orientation vectors.
            particles[i].position.resize(box.dimension);
            particles[i].orientation.resize(box.dimension);

            // Load position.
            for (unsigned int j=0;j<box.dimension;j++)
                dataFile >> particles[i].position[j];

            // Load orientation.
            // Becuase we store the quaternion elements (sin(t/2), cos(t/2)) in the
            // trajectory, we need to convert them to orientation vectors. The conversion is
            //    sin(t) = 2*sin(t/2)*cos(t/2)
            //    cos(t) = cos(t/2)^2 - sin(t/2)^2
            // where orientation is just (sin(t), cos(t))
            // We make sure that orientation is normal at the end.
            for (unsigned int j=0;j<box.dimension;j++)
                dataFile >> dummy[j];
            particles[i].orientation[0] = dummy[1]*dummy[1] - dummy[0]*dummy[0];
            particles[i].orientation[1] = 2 * dummy[0]*dummy[1];
            double norm = sqrt(particles[i].orientation[0]*particles[i].orientation[0]
                             + particles[i].orientation[1]*particles[i].orientation[1]);
            particles[i].orientation[0] /= norm;
            particles[i].orientation[1] /= norm;

            // Enforce periodic boundary conditions.
            box.periodicBoundaries(particles[i].position);

            // Calculate the particle's cell index.
            particles[i].cell = cells.getCell(particles[i]);

            // Update cell list.
            cells.initCell(particles[i].cell, particles[i]);
        }
    }
    else
    {
        std::cerr << "[ERROR] loadInitialConfiguration(): Invalid init_conf file ('"<< fileName <<"')!\n";
        exit(EXIT_FAILURE);
    }
    // Close file stream.
    dataFile.close();
}

// Maybe I should give as input TOP
void loadInitialConfigurationPatches(std::string fileName, std::vector<Particle>& particles)
{
    std::ifstream dataFile;
    std::vector<double> dummy (dimension);

    dataFile.open(fileName.c_str());
    // to read something dataFile >>
    // Check that the file is valid.
    if (dataFile.good())
    {
        for (unsigned int i=0;i<nParticles;i++)
        {
            for (unsigned int jj=0;jj<top.nPatches[particles[i].type] ;jj++)

            // Set patchstatus TO CHECK
            dataFile >> particles[i].patchstates[jj];
            // Note that as long as the file has >= patchstatus reported the reading will go well and simulation start. be careful!
        }
    }
    else
    {
     	std::cerr << "[ERROR] loadInitialConfigurationPatches(): Invalid init_conf file ('"<< fileName <<"')!\n";
        exit(EXIT_FAILURE);
    }
    // Close file stream.
    dataFile.close();
}


void parseInitialisationBlock()
{
    std::string name="initialisation";
    //Zzstd::cout << name;
    // figure out initialisation mode

    if (inp[name].contains("mode")) init_mode.assign(inp[name]["mode"]);
    //init_mode::assign()
    if (inp[name].contains(init_mode))
    {
        if (init_mode == "from_random_conf")
        {
            nTypes =inp[name][init_mode]["types"];
            top.setSize(nTypes);
            for (unsigned int i=0;i<nTypes;i++){
                unsigned int t1 = inp[name][init_mode]["particle_numbers"][i]["type"];
                top.Ni[t1] = inp[name][init_mode]["particle_numbers"][i]["N"];
                nParticles += top.Ni[t1];
                //std::cout<< nParticles << std::endl;
            }
            
            for (unsigned int k=0; k<dimension; k++)
            {
                boxSizeVec[k] = inp[name][init_mode]["box"][k];
            }
            surface_interaction=inp[name][init_mode]["surface"];
        }
        else if (init_mode == "from_init_conf")
        {
            nTypes =inp[name][init_mode]["types"];
            top.setSize(nTypes);
            for (unsigned int i=0;i<nTypes;i++){
                unsigned int t1 = inp[name][init_mode]["particle_numbers"][i]["type"];
                top.Ni[t1] = inp[name][init_mode]["particle_numbers"][i]["N"];
                nParticles += top.Ni[t1];
                //std::cout<< nParticles << std::endl;
            }
            
            for (unsigned int k=0; k<dimension; k++)
            {
             	boxSizeVec[k] = inp[name][init_mode]["box"][k];
            }
            surface_interaction=inp[name][init_mode]["surface"];
            init_conf.assign(inp[name][init_mode]["init_conf"]);
            init_patch_conf.assign(inp[name][init_mode]["init_patch_conf"]);
            restart_step_counter = inp[name][init_mode]["restart_step_counter"];
        }
    }
    else
    {
        std::cerr << "[ERROR] Initialisation mode is not defined!\n";
        exit(EXIT_FAILURE);
    }
}

void parseTopologyBlock()
{
    // Initialise topology
    // patches
    unsigned int t1, t2;
    std::string name="topology";
    for (unsigned int i=0;i<nTypes;i++){
        t1 = inp[name]["patches"][i]["type"];
        top.nPatches[t1] = inp[name]["patches"][i]["nPatches"];
        //std::cout << inp["patches"][i]["angles"] << std::endl;
        for (unsigned int j=0;j<top.nPatches[t1];j++){
            double dummy = inp[name]["patches"][i]["angles"][j];
            top.patchAngles[t1].push_back(M_PI * (dummy / 180.0));
        }
        top.patchstates[t1] = inp[name]["patches"][i]["nStates"];
        double dummy=0;
        for (unsigned int j=0;j<top.patchstates[t1];j++){
            double dummy2 = inp[name]["patches"][i]["probability"][j];
            dummy+=dummy2;
            top.registerstatus[t1].push_back(dummy);
        }
        // force a limited amount of states by putting the last element of registerstatus to a number equal to 1
        top.registerstatus[t1].push_back(1.0);
    }
    //init pair_coeff
    unsigned int counter = 0;
    for (unsigned int i=0;i<nTypes;i++){
        for (unsigned int j=0;j<nTypes;j++){
            t1 = inp[name]["pair_coeff"][counter]["type1"];
            t2 = inp[name]["pair_coeff"][counter]["type2"];
            top.epsilon[t1][t2] = inp[name]["pair_coeff"][counter]["epsilon"];
            top.delta[t1][t2] = inp[name]["pair_coeff"][counter]["delta"];
            top.sigma[t1][t2] = inp[name]["pair_coeff"][counter]["sigma"];
            top.sigma_p[t1][t2] = inp[name]["pair_coeff"][counter]["sigma_p"]; // Sigma p in case of the KFOpenSurfPatchyDisc it is the maximum angle for the interaction
            top.rcut[t1][t2] = inp[name]["pair_coeff"][counter]["rcut"];
            top.shift[t1][t2] = inp[name]["pair_coeff"][counter]["shift"];
            counter += 1;
            if (interactionEnergy<top.epsilon[t1][t2]) interactionEnergy=top.epsilon[t1][t2];
            if (interactionRange<top.rcut[t1][t2]) interactionRange=top.rcut[t1][t2];
        }
    }
}

void parseSimulationParamBlock()
{
    std::string name = "simulation parameters";
    if (inp[name].contains("output_every")) output_every = inp[name]["output_every"];
    if (inp[name].contains("last_conf")) last_conf.assign(inp[name]["last_conf"]);
    if (inp[name].contains("trajectory")) trajectory.assign(inp[name]["trajectory"]);
    if (inp[name].contains("log_file")) log_file.assign(inp[name]["log_file"]);
    if (inp[name].contains("seed")) seed = inp[name]["seed"];
    if (inp[name].contains("interaction")) interaction.assign(inp[name]["interaction"]);
    sweeps = inp[name]["sweeps"];
}


void parseInputFile()
{
    parseInitialisationBlock();
    parseTopologyBlock();
    parseSimulationParamBlock();
}

#endif
