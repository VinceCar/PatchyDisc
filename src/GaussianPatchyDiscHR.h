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

  Adapted from PatchyDisc.h
*/

#ifndef _GAUSSIANPATCHYDISCHR_H
#define _GAUSSIANPATCHYDISCHR_H

#include "Model.h"

/*! \file GaussianPatchyDiscHR.h
*/

// FORWARD DECLARATIONS
class Top;

//! Class defining the Patchy-Disc potential.
class GaussianPatchyDiscHR : public Model
{
public:
    //! Constructor.
    /*! \param box_
            A reference to the simulation box object.

        \param particles_
            A reference to the particle list.

        \param cells_
            A reference to the cell list object.

        \param maxInteractions_
            The maximum number of interactions per particle (number of patches).

        \param interactionEnergy_
            The potential energy scale (in units of kBT).

        \param interactionRange_
            The potential cut-off distance (patch diameter).
     */
    GaussianPatchyDiscHR(Box&, std::vector<Particle>&, CellList&, Top&, unsigned int, double, double);
    //! Calculate the pair energy between two particles.
    /*! \param particle1
            The index of the first particle.

        \param position1
            The position vector of the first particle.

        \param orientation1
            The orientation vector of the first particle.

        \param particle2
            The index of the second particle.

        \param position2
            The position vector of the second particle.

        \param orientation2
            The orientation vector of the second particle.

        \return
            The pair energy between particles 1 and 2.
     */
    double computePairEnergy(unsigned int, const double*, const double*, unsigned int, const double*, const double*);

    //! Determine the interactions for a given particle.
    /*! \param particle
            The particle index.

        \param position
            The position vector of the particle.

        \param orientation
            The orientation vector of the particle.

        \param interactions
            An array to store the indices of neighbours with which the particle interacts.

        \return
            The number of interactions.
     */
    unsigned int computeInteractions(unsigned int, const double*, const double*, unsigned int*);

    Top& top;                           //!< A reference to the Top.

private:
    std::vector<unsigned int> idx2type;
    std::vector<std::vector<double>> rcut_sq;
    std::vector<std::vector<double>> lj_shift;
    std::vector<std::vector<double>> twoSigmapSq;
    std::vector<std::vector<double>> sigmaSq;

    // std::vector<double> expTable;
    // double dx = 5e-4;
    // double max_x = 10;
    // unsigned int tableSize = max_x / dx + 1.0;
};

#endif  /* _GAUSSIANPATCHYDISCHR_H */