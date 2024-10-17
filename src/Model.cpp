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

#include <cstdlib>
#include <iostream>
#include <limits>
#include <cmath>
#include <csignal>
#include "Box.h"
#include "CellList.h"
#include "Model.h"
#include "Particle.h"
#include "Top.h"
double INF = std::numeric_limits<double>::infinity();

Model::Model(
    Box& box_,
    std::vector<Particle>& particles_,
    CellList& cells_,
    Top& top_,
    std::vector<unsigned int>& bonds_,
    unsigned int maxInteractions_,
    double interactionEnergy_,
    double interactionRange_) :

    box(box_),
    particles(particles_),
    cells(cells_),
    top(top_),
    bonds(bonds_),
    maxInteractions(maxInteractions_),
    interactionEnergy(interactionEnergy_),
    interactionRange(interactionRange_)
{
    // Work out squared cut-off distance.
    squaredCutOffDistance = interactionRange * interactionRange;
}

#ifndef ISOTROPIC
double Model::computeEnergy(unsigned int particle, unsigned int cell_particle, const double* position, const double* orientation, const unsigned int* status, std::vector<int>& id_contacts, std::vector<unsigned int>& bonds)
#else
double Model::computeEnergy(unsigned int particle, unsigned int cell_particle, const double* position, std::vector<int>& id_contacts, std::vector<unsigned int>& bonds)
#endif
{
    // N.B. This method is somewhat redundant since the same functionality
    // could be achieved by using a combination of the computeInteractions
    // and model specific computePairEnergy methods.

    // Energy counter.
    double energy=0;
    double energy_single=0;
    bonds[particle]=0;
    // Check all neighbouring cells including same cell.
    bool flag_inf=false;
    for (unsigned int i=0;i<cells.getNeighbours();i++)
    {
        // Cell index.
        unsigned int cell = cells[cell_particle].neighbours[i];

        // Check all particles within cell.
        for (unsigned int j=0;j<cells[cell].tally;j++)
        {
            // Index of neighbouring particle.
            unsigned int neighbour = cells[cell].particles[j];

            // Make sure the particles are different.
            if (neighbour != particle)
            {
                // Calculate model specific pair energy.
#ifndef ISOTROPIC
                energy_single = computePairEnergy(particle, position, orientation, status,
                          neighbour, &particles[neighbour].position[0],
                          &particles[neighbour].orientation[0], &particles[neighbour].patchstates[0]);
#else
                energy_single = computePairEnergy(particle, position,
                          neighbour, &particles[neighbour].position[0]);
#endif
                if (energy_single<0)
                {
                   bonds[particle]++;
                   id_contacts.push_back(neighbour);
                }
                if (energy_single>=1e8)
                {
                    flag_inf=true;
                }
                energy+=energy_single;
            }
        }
    }
    if (flag_inf) {
       energy=INF;
    }
    return energy;
}

#ifndef ISOTROPIC
double Model::computePairEnergy(unsigned int particle1, const double* position1, const double* orientation1, const unsigned int* status1,
    unsigned int particle2, const double* position2, const double* orientation2, const unsigned int* status2)
#else
double Model::computePairEnergy(unsigned int particle1,
    const double* position1, unsigned int particle2, const double* position2)
#endif
{
    std::cerr << "[ERROR] Model: Virtual function Model::computePairEnergy() must be defined.\n";
    exit(EXIT_FAILURE);
}

#ifndef ISOTROPIC
unsigned int Model::computeInteractions(unsigned int particle,
    const double* position, const double* orientation, unsigned int* interactions)
#else
unsigned int Model::computeInteractions(unsigned int particle,
    const double* position, unsigned int* interactions)
#endif
{
    // Interaction counter.
    unsigned int nInteractions = 0;

    // Check all neighbouring cells including same cell.
    for (unsigned int i=0;i<cells.getNeighbours();i++)
    {
        // Cell index.
        unsigned int cell = cells[particles[particle].cell].neighbours[i];

        // Check all particles within cell.
        for (unsigned int j=0;j<cells[cell].tally;j++)
        {
            // Index of neighbouring particle.
            unsigned int neighbour = cells[cell].particles[j];

            // Make sure the particles are different.
            if (neighbour != particle)
            {
                std::vector<double> sep(box.dimension);

                // Compute separation.
                for (unsigned int k=0;k<box.dimension;k++)
                    sep[k] = position[k] - particles[neighbour].position[k];

                // Enforce minimum image.
                box.minimumImage(sep);

                double normSqd = 0;

                // Calculate squared norm of vector.
                for (unsigned int k=0;k<box.dimension;k++)
                    normSqd += sep[k]*sep[k];

                // Particles interact.
                if (normSqd < squaredCutOffDistance)
                {
                    if (nInteractions == maxInteractions)
                    {
                        std::cerr << "[ERROR] Model: Maximum number of interactions exceeded!\n";
                        exit(EXIT_FAILURE);
                    }

                    interactions[nInteractions] = neighbour;
                    nInteractions++;
                }
            }
        }
    }

    return nInteractions;
}

#ifndef ISOTROPIC
void Model::applyPostMoveUpdates(unsigned int particle, const double* position, const double* orientation, const unsigned int* status)
#else
void Model::applyPostMoveUpdates(unsigned int particle, const double* position)
#endif
{
    // Copy coordinates/orientations.
    for (unsigned int i=0;i<box.dimension;i++)
    {
        particles[particle].position[i] = position[i];
#ifndef ISOTROPIC
        {
        particles[particle].orientation[i] = orientation[i];
        particles[particle].patchstates[i] = status[i];
        }
#endif
    }

    // Calculate the particle's cell index.
    unsigned int newCell = cells.getCell(particles[particle]);

    // Update cell lists if necessary.
    if (particles[particle].cell != newCell)
        cells.updateCell(newCell, particles[particle], particles);
}

double Model::getEnergy()
{
    double energy = 0;
    std::vector<int> id_contacts;

    for (unsigned int i=0;i<particles.size();i++)
#ifndef ISOTROPIC
        energy += computeEnergy(i, particles[i].cell, &particles[i].position[0], &particles[i].orientation[0], &particles[i].patchstates[0], id_contacts, bonds);
#else
        energy += computeEnergy(i, particles[i].cell, &particles[i].position[0], id_contacts, bonds);
#endif

    return energy/(2*particles.size());
}
