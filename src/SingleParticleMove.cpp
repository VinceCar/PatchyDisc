/*
  Copyright (c) 2015-2016 Lester Hedges <lester.hedges+vmmc@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/
#include <iostream>
#include <csignal>
#include "Box.h"
#include "CellList.h"
#include "Model.h"
#include "Particle.h"
#include "SingleParticleMove.h"
#include "Top.h"

SingleParticleMove::SingleParticleMove(
    MersenneTwister& rng_,
    Model* model_,
    double surf_interaction_,
    double maxTrialTranslation_,
    double maxTrialRotation_,
    double probTranslate_,
    double probChange_,
    bool isIsotropic_) :

    rng(rng_),
    model(model_),
    surf_interaction(surf_interaction_),
    maxTrialTranslation(maxTrialTranslation_),
    maxTrialRotation(maxTrialRotation_),
    probTranslate(probTranslate_),
    probChange(probChange_),
    isIsotropic(isIsotropic_)

{
    // Check dimensionality.
    if (model->box.dimension == 3) is3D = true;
    else is3D = false;

    // Allocate memory.
    moveParams.trialVector.resize(model->box.dimension);

    // Ignore rotations if potential is isotropic.
    if (isIsotropic) probTranslate = 1.0;
    std::cout << "# Initialised MC" << std::endl;

}

void SingleParticleMove::step(const int nSteps)
{
    for (int i=0;i<nSteps;i++)
        step();
}

void SingleParticleMove::operator ++ (const int)
{
    step();
}

void SingleParticleMove::operator += (const int nSteps)
{
    step(nSteps);
}

void SingleParticleMove::step()
{
    // Increment number of attempted moves.
    nAttempts++;

    // Propose a trial move.
    if (nAttempts-2*nAccepts<=0)
    {
    	proposeMove();
    	// Check whether move was accepted.
    	if (accept())
    	{
            // Increment number of accepted moves.
            nAccepts++;
            if (moveParams.isChange)
            {
                nChanges++;
            }
            else if (!moveParams.isRotation)
            {
                // Increment number of rotations.
                nRotations += moveParams.isRotation;
                unsigned int oldCell = moveParams.preMoveParticle.cell;
                unsigned int newCell = model->particles[moveParams.seed].cell;
                // update cell list
                if (oldCell != newCell)
       	        {
                    model->particles[moveParams.seed].cell = oldCell;
                    model->cells.updateCell(newCell, model->particles[moveParams.seed], model->particles);
            	}
            } else {
                // Increment number of rotations.
                nRotations += moveParams.isRotation;
            }
    	} else {
        		// Revert particle to pre-move state.
            model->particles[moveParams.seed] = moveParams.preMoveParticle;
            model->bonds = moveParams.preMoveBonds;
    	}
    } else {
        proposeMove_post();
        // Check whether move was accepted.
        if (accept())
        {
            // Increment number of accepted moves.
            nAccepts++;
            unsigned int oldCell = model->particles[moveParams.seed].cell;
            model->particles[moveParams.seed] = moveParams.postMoveParticle;
            model->bonds = moveParams.postMoveBonds;

            if (moveParams.isChange)
            {
                nChanges++;
            } else if (!moveParams.isRotation) {
                // Increment number of rotations.
                nRotations += moveParams.isRotation;
                unsigned int newCell = moveParams.postMoveParticle.cell;
                // update cell list
                if (oldCell != newCell)
                {
                    model->particles[moveParams.seed].cell = oldCell;
                    model->cells.updateCell(newCell, model->particles[moveParams.seed], model->particles);
                }
            } else {
                // Increment number of rotations.
                nRotations += moveParams.isRotation;
            }
        }

    }
}

unsigned long long SingleParticleMove::getAttempts() const
{
    return nAttempts;
}

unsigned long long SingleParticleMove::getAccepts() const
{
    return nAccepts;
}

unsigned long long SingleParticleMove::getRotations() const
{
    return nRotations;
}

unsigned long long SingleParticleMove::getChanges() const
{
    return nChanges;
}
void SingleParticleMove::reset()
{
    nAttempts = nAccepts = nRotations = nChanges = 0;
}

void SingleParticleMove::proposeMove()
{
    // Variables for the surface potentials, they will be moved in another section
    std::vector<bool> init_contacts(model->particles.size(),false); // Flag vector of the initial contacts
    std::vector<bool> final_contacts(model->particles.size(),false); // Flag vector of the  final contacts
    std::vector<int> init_id_contacts; // ids for the initial contacts
    std::vector<int> final_id_contacts; //ids for the final contacts
    final_id_contacts.reserve(model->maxInteractions);
    init_id_contacts.reserve(model->maxInteractions);
    int delta_contacts=0;
   // Choose a seed particle.
    moveParams.seed = rng.integer(0, model->particles.size()-1);
    double changenergy=0;
    energyChange=0;
    // Choose a random point on the surface of the unit sphere/circle.
    for (unsigned int i=0;i<model->box.dimension;i++)
        moveParams.trialVector[i] = rng.normal();

    // Calculate pre-move energy.
#ifndef ISOTROPIC
    initialEnergy = model->computeEnergy(moveParams.seed, model->particles[moveParams.seed].cell,
        &model->particles[moveParams.seed].position[0],
        &model->particles[moveParams.seed].orientation[0],
        &model->particles[moveParams.seed].patchstates[0], init_id_contacts, model->bonds);
    double initial_state=0;
    for (unsigned int ii=0;ii< model->particles[moveParams.seed].patchstates.size();ii++) {
        initial_state+=(model->particles[moveParams.seed].patchstates[ii]>0);
    }
#else
    initialEnergy = model->computeEnergy(moveParams.seed, model->particles[moveParams.seed].cell, 
        &model->particles[moveParams.seed].position[0], init_id_contacts, model->bonds);
    double initial_state=0;
#endif

    // Assign a flag there
    // Normalise the trial vector.
    double norm = computeNorm(moveParams.trialVector);
    for (unsigned int i=0;i<model->box.dimension;i++)
        moveParams.trialVector[i] /= norm;
     moveParams.isChange = false;
     moveParams.isRotation = false;


    // Choose the move type.
    if (rng() < probChange)
    {
        // Change patch type
        moveParams.isChange = true;
        moveParams.stepSize = rng();
    } else {
        if (rng() < probTranslate)
        {
            // Translation.
            // Scale step-size to uniformly sample unit sphere/circle.
            if (is3D) moveParams.stepSize = maxTrialTranslation*std::pow(rng(), 1.0/3.0);
            else moveParams.stepSize = maxTrialTranslation*std::pow(rng(), 1.0/2.0);
        } else {
            // Rotation.
            moveParams.isRotation = true;
            moveParams.stepSize = maxTrialRotation*(2.0*rng()-1.0);
        }
    }

    // Store initial coordinates/orientation.
    moveParams.preMoveParticle = model->particles[moveParams.seed];
    moveParams.preMoveBonds=model->bonds;
    // Execute the move.
    if (moveParams.isChange)
    {
       changenergy = changepatch(model->particles[moveParams.seed].patchstates, model->top.registerstatus[model->particles[moveParams.seed].type], moveParams.stepSize);
    } else if (!moveParams.isRotation) // Translation.
    {
        for (unsigned int i=0;i<model->box.dimension;i++)
        {
            model->particles[moveParams.seed].position[i] += moveParams.stepSize*moveParams.trialVector[i];

        // Apply periodic boundary conditions.
        model->box.periodicBoundaries(model->particles[moveParams.seed].position);

        // Work out new cell index.
        model->particles[moveParams.seed].cell = model->cells.getCell(model->particles[moveParams.seed]);
        }
    } else {
        std::vector<double> vec(model->box.dimension);

        // Calculate orientation rotation vector.
        if (is3D) rotate3D(model->particles[moveParams.seed].orientation, moveParams.trialVector, vec, moveParams.stepSize);
        else rotate2D(model->particles[moveParams.seed].orientation, vec, moveParams.stepSize);

        // Update orientation.
        for (unsigned int i=0;i<model->box.dimension;i++)
        {
            model->particles[moveParams.seed].orientation[i] += vec[i];
        }
    }

    // Calculate post-move energy.
#ifndef ISOTROPIC
    double finalEnergy = model->computeEnergy(moveParams.seed, model->particles[moveParams.seed].cell, 
        &model->particles[moveParams.seed].position[0],
        &model->particles[moveParams.seed].orientation[0],
        &model->particles[moveParams.seed].patchstates[0], final_id_contacts, model->bonds);
    double final_state=0;
    for (unsigned int ii=0;ii< model->particles[moveParams.seed].patchstates.size();ii++) {
        final_state+=(model->particles[moveParams.seed].patchstates[ii]>0);
    }

#else
    double finalEnergy = model->computeEnergy(moveParams.seed, model->particles[moveParams.seed].cell, 
        &model->particles[moveParams.seed].position[0], final_id_contacts, model->bonds);
    double final_state=0;
#endif
   // Setup the variables for surface interaction
    for (int index : init_id_contacts) {
        if (index < init_contacts.size()) {
            init_contacts[index] = true;
        } else {
            std::cout << "Error, out of boundary access when setting the neighbours" << std::endl;
        }
    }
    for (int index : final_id_contacts) {
        if (index < final_contacts.size()) {
            final_contacts[index] = true;
        } else {
            std::cout << "Error, out of boundary access when setting the neighbours" << std::endl;
        }
    }

    for (int index : final_id_contacts) {
        delta_contacts = final_contacts[index] - init_contacts[index];
        if (delta_contacts > 0) {
            if (model->bonds[index] == 0) {
                energyChange -= surf_interaction;
            }
            model->bonds[index]++;
        }
    }
    for (int index : init_id_contacts) {
        delta_contacts = final_contacts[index] - init_contacts[index];
        if (delta_contacts < 0) {
            model->bonds[index]--;
            if (model->bonds[index] == 0) {
                energyChange += surf_interaction;
            }
	}
    }

   // calculate the total energy change
    energyChange += (finalEnergy - initialEnergy) + changenergy ;
}


void SingleParticleMove::proposeMove_post()
{
    // Variables for the surface potentials, they will be moved in another section
    std::vector<bool> init_contacts(model->particles.size(),false); // Flag vector of the initial contacts
    std::vector<bool> final_contacts(model->particles.size(),false); // Flag vector of the  final contacts
    std::vector<int> init_id_contacts; // ids for the initial contacts
    std::vector<int> final_id_contacts; //ids for the final contacts
    final_id_contacts.reserve(model->maxInteractions);
    init_id_contacts.reserve(model->maxInteractions);
    energyChange=0;
    int delta_contacts=0;
   // Choose a seed particle.
    moveParams.seed = rng.integer(0, model->particles.size()-1);
    double changenergy=0;
    // Choose a random point on the surface of the unit sphere/circle.
    for (unsigned int i=0;i<model->box.dimension;i++)
        moveParams.trialVector[i] = rng.normal();

    // Calculate pre-move energy.
#ifndef ISOTROPIC
    initialEnergy = model->computeEnergy(moveParams.seed, model->particles[moveParams.seed].cell, 
        &model->particles[moveParams.seed].position[0],
        &model->particles[moveParams.seed].orientation[0],
        &model->particles[moveParams.seed].patchstates[0], init_id_contacts, model->bonds);
    double initial_state=0;
    for (unsigned int ii=0;ii< model->particles[moveParams.seed].patchstates.size();ii++) {
        initial_state+=(model->particles[moveParams.seed].patchstates[ii]>0);
    }
#else
    initialEnergy = model->computeEnergy(moveParams.seed, model->particles[moveParams.seed].cell,
        &model->particles[moveParams.seed].position[0], init_id_contacts, model->bonds);
    double initial_state=0;
#endif

    // Assign a flag there
    // Normalise the trial vector.
    double norm = computeNorm(moveParams.trialVector);
    for (unsigned int i=0;i<model->box.dimension;i++)
        moveParams.trialVector[i] /= norm;
     moveParams.isChange = false;
     moveParams.isRotation = false;


    // Choose the move type.
    if (rng() < probChange)
    {
     	// Change patch type
        moveParams.isChange = true;
        moveParams.stepSize = rng();
    } else {
	if (rng() < probTranslate)
        {
            // Translation.
            // Scale step-size to uniformly sample unit sphere/circle.
            if (is3D) moveParams.stepSize = maxTrialTranslation*std::pow(rng(), 1.0/3.0);
            else moveParams.stepSize = maxTrialTranslation*std::pow(rng(), 1.0/2.0);
        } else {
            // Rotation.
            moveParams.isRotation = true;
            moveParams.stepSize = maxTrialRotation*(2.0*rng()-1.0);
        }
    }

    // Store initial coordinates/orientation.
    moveParams.postMoveParticle = model->particles[moveParams.seed];
    moveParams.postMoveBonds=model->bonds;
    // Execute the move.
    if (moveParams.isChange)
    {
       changenergy = changepatch(moveParams.postMoveParticle.patchstates, model->top.registerstatus[model->particles[moveParams.seed].type], moveParams.stepSize);
    } else if (!moveParams.isRotation) // Translation.
    {
     	for (unsigned int i=0;i<model->box.dimension;i++)
        {
            moveParams.postMoveParticle.position[i] += moveParams.stepSize*moveParams.trialVector[i];

        // Apply periodic boundary conditions.
        model->box.periodicBoundaries(moveParams.postMoveParticle.position);

        // Work out new cell index.
        moveParams.postMoveParticle.cell = model->cells.getCell(moveParams.postMoveParticle);
        }
    } else {
	std::vector<double> vec(model->box.dimension);

        // Calculate orientation rotation vector.
        if (is3D) rotate3D(moveParams.postMoveParticle.orientation, moveParams.trialVector, vec, moveParams.stepSize);
        else rotate2D(moveParams.postMoveParticle.orientation, vec, moveParams.stepSize);

        // Update orientation.
        for (unsigned int i=0;i<model->box.dimension;i++)
        {
            moveParams.postMoveParticle.orientation[i] += vec[i];
        }
    }

    // Calculate post-move energy.
#ifndef ISOTROPIC
    double finalEnergy = model->computeEnergy(moveParams.seed, moveParams.postMoveParticle.cell,   
        &moveParams.postMoveParticle.position[0],
        &moveParams.postMoveParticle.orientation[0],
        &moveParams.postMoveParticle.patchstates[0], final_id_contacts, moveParams.postMoveBonds);
    double final_state=0;
    for (unsigned int ii=0;ii< moveParams.postMoveParticle.patchstates.size();ii++) {
        final_state+=(moveParams.postMoveParticle.patchstates[ii]>0);
    }

#else
    double finalEnergy = model->computeEnergy(moveParams.seed, moveParams.postMoveParticle.cell, 
        &moveParams.postMoveParticle.position[0], final_id_contacts, moveParams.postMoveBonds);
    double final_state=0;
#endif

   // Setup the variables for surface interaction
    for (int index : init_id_contacts) {
        if (index < init_contacts.size()) {
            init_contacts[index] = true;
        } else {
            std::cout << "Error, out of boundary access when setting the neighbours" << std::endl;
        }
    }
    for (int index : final_id_contacts) {
        if (index < final_contacts.size()) {
            final_contacts[index] = true;
        } else {
            std::cout << "Error, out of boundary access when setting the neighbours" << std::endl;
        }
    }

    for (int index : final_id_contacts) {
        delta_contacts = final_contacts[index] - init_contacts[index];
        if (delta_contacts > 0) {
            if (model->bonds[index] == 0) {
                energyChange -= surf_interaction;
            }
            model->bonds[index]++;
        }
    }
    for (int index : init_id_contacts) {
        delta_contacts = final_contacts[index] - init_contacts[index];
        if (delta_contacts < 0) {
            model->bonds[index]--;
            if (model->bonds[index] == 0) {
                energyChange += surf_interaction;
            }
        }
    }

   // calculate the total energy change
    energyChange += (finalEnergy - initialEnergy) + changenergy ;
}

bool SingleParticleMove::accept()
{
    if (energyChange == 0) return true;
    if (energyChange == INF) return false;
    if ((rng() < exp(-energyChange))) return true;
    else return false;
}

void SingleParticleMove::rotate3D(std::vector<double>& v1, std::vector<double>& v2, std::vector<double>& v3, double angle)
{
    double c = cos(angle);
    double s = sin(angle);

    double v1Dotv2 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];

    v3[0] = ((v1[0] - v2[0]*v1Dotv2))*(c - 1) + (v2[2]*v1[1] - v2[1]*v1[2])*s;
    v3[1] = ((v1[1] - v2[1]*v1Dotv2))*(c - 1) + (v2[0]*v1[2] - v2[2]*v1[0])*s;
    v3[2] = ((v1[2] - v2[2]*v1Dotv2))*(c - 1) + (v2[1]*v1[0] - v2[0]*v1[1])*s;
}

void SingleParticleMove::rotate2D(std::vector<double>& v1, std::vector<double>& v2, double angle)
{
    double c = cos(angle);
    double s = sin(angle);

    v2[0] = (v1[0]*c - v1[1]*s) - v1[0];
    v2[1] = (v1[0]*s + v1[1]*c) - v1[1];
}

double SingleParticleMove::changepatch(std::vector<unsigned int>& v1, std::vector<double>& v2, double probability)
{
    int numberofpatches = v1.size();
    int patchtochange = (int)round(numberofpatches*rng()-0.5);
    unsigned int ii = (unsigned int) round((v2.size()-1.5)*probability);
    double previous=v2[v1[patchtochange]];
    if (v1[patchtochange]>0)
    {
           previous=previous-v2[v1[patchtochange]-1];
    }
    if (ii>=v1[patchtochange])
    {
        ii++;
    }
    v1[patchtochange] = ii;
    if (v1[patchtochange]>0)
    {
         previous=previous/(v2[v1[patchtochange]]-v2[v1[patchtochange]-1]);
    } else {
         previous=previous/(v2[v1[patchtochange]]);
    }
    return log(previous);
}

double SingleParticleMove::computeNorm(std::vector<double>& vec)
{
    double normSquared = 0;

    for (unsigned int i=0;i<vec.size();i++)
        normSquared += vec[i]*vec[i];

    return sqrt(normSquared);
}
