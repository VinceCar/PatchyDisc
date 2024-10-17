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

  Adapted from PatchyDisc.cpp
*/

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <csignal>
#include <iomanip>      // std::setw
#include <limits>
#include <limits.h>     // CHAR_BIT
#include <sstream>
#include <stdint.h>     // uint64_t
#include "Box.h"
#include "CellList.h"
#include "Particle.h"
#include "Top.h"
#include "KFOpenSurfPatchyDisc.h"
#include "Utils.h"

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif

KFOpenSurfPatchyDisc::KFOpenSurfPatchyDisc(
    Box& box_,
    std::vector<Particle>& particles_,
    CellList& cells_,
    Top& top_,
    std::vector<unsigned int>& bonds_,
    unsigned int maxInteractions_,
    double interactionEnergy_,
    double interactionRange_) :
    Model(box_, particles_, cells_, top_, bonds_, maxInteractions_, interactionEnergy_, interactionRange_),
    top(top_)
{
#ifdef ISOTROPIC
    std::cerr << "[ERROR] KFOpenSurfPatchyDisc: Cannot be used with isotropic VMMC library!\n";
    exit(EXIT_FAILURE);
#endif

    // Check dimensionality.
    if (box.dimension != 2)
    {
        std::cerr << "[ERROR] KFOpenSurfPatchyDisc: Model only valid in two dimensions!\n";
        exit(EXIT_FAILURE);
    }
    //init all variables and look-up tables
    for (unsigned int i=0;i<particles.size();i++)
    {
        idx2type.push_back(particles[i].type);
    }
    rcut_sq.resize(top.nTypes, std::vector<double>(top.nTypes));
    lj_shift.resize(top.nTypes, std::vector<double>(top.nTypes));
    twoSigmapSq.resize(top.nTypes, std::vector<double>(top.nTypes)); // Note it is not twoSigmapsquare anymore, it is just the cosine of the maximum angle
    sigmaSq.resize(top.nTypes, std::vector<double>(top.nTypes));
    for (unsigned int i=0;i<top.nTypes;i++)
    {
        for (unsigned int j=0;j<top.nTypes;j++)
        {
            rcut_sq[i][j] = top.rcut[i][j] * top.rcut[i][j];
            lj_shift[i][j] = pow(top.sigma[i][j]/top.rcut[i][j], 12) - pow(top.sigma[i][j]/top.rcut[i][j], 6);
            lj_shift[i][j] *= -4 * top.epsilon[i][j] * top.shift[i][j];
            twoSigmapSq[i][j] = cos(top.sigma_p[i][j]);
            sigmaSq[i][j] = top.sigma[i][j] * top.sigma[i][j];
         }
    }

    // expTable.resize(tableSize);
    // for (int i=0; i<tableSize; i++)
    // {
    //     expTable[i] = exp(-i*dx);
    // }
}

double KFOpenSurfPatchyDisc::computePairEnergy(unsigned int particle1, const double* position1,
    const double* orientation1, const unsigned int* status1, unsigned int particle2, const double* position2, const double* orientation2, const unsigned int* status2)
{

    unsigned int t1 = idx2type[particle1];
    unsigned int t2 = idx2type[particle2];
    // Separation vector.
    std::vector<double> sep(2);
    // Calculate disc separation.
    sep[0] = position2[0] - position1[0];
    sep[1] = position2[1] - position1[1];
    // Enforce minimum image.
    box.minimumImage(sep);

    // Calculate squared norm of vector.
    double normSqd = sep[0]*sep[0] + sep[1]*sep[1];
    // Particles interact.
    double energyLJ = 0;
    if (normSqd < sigmaSq[t1][t2])
    {
        energyLJ = 1e8;
        return energyLJ;

    }
    else if (normSqd < rcut_sq[t1][t2])
    {
      energyLJ = -top.epsilon[t1][t2];
    }
    else 
    {
        return energyLJ; 
    }
    // Calculate the angular modulation of the KF interaction.
    // examine all patch pairs and take the one that have minimum distance.
    double modulation, p1Angle, p2Angle;
    double cosalpha=1;
    double sinalpha=0;
    double distance_max=0;
    double distance=0;
    double positionpatch1=0;
    double positionpatch2=0;
    double angletouse=2*M_PI;
    double scalarproduct=0;
    double norm1=0;
    double midpointposition1=0;
    double midpointposition2=0;
    modulation = 0;
    p1Angle = 0;
    p2Angle = 0;

    //midpointpoisition
    midpointposition1=sep[0]/2;
    midpointposition2=sep[1]/2;
    // finding the closest patch
    distance_max=2*sqrt(normSqd);
    distance=distance_max;
    double modulation1=1;
    for (unsigned int i=0;i<top.nPatches[t1];i++)
    {
        cosalpha=orientation1[0];
        sinalpha=orientation1[1];
        p1Angle = top.patchAngles[t1][i];
        positionpatch1=cosalpha*cos(p1Angle)-sinalpha*sin(p1Angle);
        positionpatch1=positionpatch1*sqrt(sigmaSq[t1][t2])/2;
        positionpatch2=cosalpha*sin(p1Angle)+cos(p1Angle)*sinalpha;
        positionpatch2=positionpatch2*sqrt(sigmaSq[t1][t2])/2;
        distance=sqrt((positionpatch1-midpointposition1)*(positionpatch1-midpointposition1)+(positionpatch2-midpointposition2)*(positionpatch2-midpointposition2));
        if (distance_max>distance) {
            distance_max=distance; 
            scalarproduct=positionpatch1*midpointposition1+positionpatch2*midpointposition2;
            norm1=sqrt(positionpatch1*positionpatch1+positionpatch2*positionpatch2);
            scalarproduct=2*(scalarproduct/sqrt(normSqd)/norm1); // not the scalar product anymore but lazy to name another variable it is the cosine
	    angletouse=scalarproduct;
            modulation=status1[i]; // monitor the patch state
        }
     }
     p1Angle=angletouse; // This is actually the cosine of the angle to use
     modulation1=modulation;
     distance_max=2*sqrt(normSqd);
     for (unsigned int j=0;j<top.nPatches[t2];j++)
     {
         cosalpha=orientation2[0];
         sinalpha=orientation2[1];
         p2Angle = top.patchAngles[t2][j];
         positionpatch1=cosalpha*cos(p2Angle)-sinalpha*sin(p2Angle);
         positionpatch2=cosalpha*sin(p2Angle)+cos(p2Angle)*sinalpha;
         positionpatch1=positionpatch1*sqrt(sigmaSq[t1][t2])/2;
         positionpatch2=positionpatch2*sqrt(sigmaSq[t1][t2])/2;
         positionpatch1=positionpatch1+sep[0];
         positionpatch2=positionpatch2+sep[1];
         distance=sqrt((positionpatch1-midpointposition1)*(positionpatch1-midpointposition1)+(positionpatch2-midpointposition2)*(positionpatch2-midpointposition2));
         if (distance_max>distance) {
             distance_max=distance;
             positionpatch1=positionpatch1-sep[0]; // vector from center to patch
             positionpatch2=positionpatch2-sep[1]; // vector from center to patch
             // scalar product inverting the direction
             scalarproduct=-(positionpatch1*midpointposition1+positionpatch2*midpointposition2);
             norm1=sqrt(sigmaSq[t1][t2])/2;
             scalarproduct=2*(scalarproduct/sqrt(normSqd)/norm1); // contains the cosine of the angle of interest
            angletouse=scalarproduct;
            modulation=modulation1*status2[j]; // apply the patch state
         }
      }
      p2Angle=angletouse; // This is actually the cosine of the angle
      if (p1Angle<=twoSigmapSq[t1][t2])
      {
          modulation = 0;
      }
      if (p2Angle<=twoSigmapSq[t1][t2])
      {
         modulation = 0;
      }
      return energyLJ*modulation;
}

unsigned int KFOpenSurfPatchyDisc::computeInteractions(unsigned int particle,
    const double* position, const double* orientation, unsigned int* interactions)
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

                if (normSqd < rcut_sq[0][0])
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


// unsigned int KFOpenSurfPatchyDisc::computePatchyBonds(unsigned int particle,
//     const double* position, const double* orientation, unsigned int* interactions)
// {
//     // This is similar to computeInteractions, but returns a vector with all bonded patchy particles

//     // Interaction counter.
//     unsigned int nInteractions = 0;

//     // Check all neighbouring cells including same cell.
//     for (unsigned int i=0;i<cells.getNeighbours();i++)
//     {
//         // Cell index.
//         unsigned int cell = cells[particles[particle].cell].neighbours[i];

//         // Check all particles within cell.
//         for (unsigned int j=0;j<cells[cell].tally;j++)
//         {
//             // Index of neighbouring particle.
//             unsigned int neighbour = cells[cell].particles[j];

//             // Make sure the particles are different.
//             if (neighbour != particle)
//             {
//                 std::vector<double> sep(box.dimension);

//                 // Compute separation.
//                 for (unsigned int k=0;k<box.dimension;k++)
//                     sep[k] = position[k] - particles[neighbour].position[k];

//                 // Enforce minimum image.
//                 box.minimumImage(sep);

//                 double normSqd = 0;

//                 // Calculate squared norm of vector.
//                 for (unsigned int k=0;k<box.dimension;k++)
//                     normSqd += sep[k]*sep[k];

//                 // Particles interact.

//                 if (normSqd < rcut_sq[idx2type[particle]][idx2type[neighbour]])
//                 {
//                     if (nInteractions == maxInteractions)
//                     {
//                         std::cerr << "[ERROR] Model: Maximum number of interactions exceeded!\n";
//                         exit(EXIT_FAILURE);
//                     }

//                     interactions[nInteractions] = neighbour;
//                     nInteractions++;
//                 }
//             }
//         }
//     }

//     return nInteractions;
// }

