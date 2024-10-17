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

#ifndef _MODEL_H
#define _MODEL_H

#include <vector>

/*! \file Model.h
*/

// FORWARD DECLARATIONS
class  Box;
class  CellList;
struct Particle;
// FORWARD DECLARATIONS
class Top;

// Global infinity constant for hard core repulsions.
extern double INF;

//! Base class defining general access functions for the model potential
//! and a virtual interface for model specific pair energies.
class Model
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
            The maximum number of interactions per particle.

        \param interactionEnergy_
            The interaction energy (in units of kBT).

        \param interactionRange_
            The interaction range (in units of the particle diameter).
     */
    Model(Box&, std::vector<Particle>&, CellList&, Top&, std::vector<unsigned int>&, unsigned int, double, double);

    virtual ~Model() = default;
    Box& box;                           //!< A reference to the simulation box.
    std::vector<Particle>& particles;   //!< A reference to the particle list.
    CellList& cells;                    //!< A reference to the cell list.
    Top& top;                           //!< A reference to the topology
    std::vector<unsigned int>& bonds;
    //! Calculate the total interaction energy felt by a particle.
    /*! \param index
            The particle index.

        \param position
            The position vector of the particle.

        \param orientation
            The orientation vector of the first particle.
        \param status
            The status vector of the patches
        \return
            The total interaction energy.
     */
#ifndef ISOTROPIC
    virtual double computeEnergy(unsigned int, unsigned int, const double*, const double*, const unsigned int*, std::vector<int>&,  std::vector<unsigned int>&);
#else
    virtual double computeEnergy(unsigned int, unsigned int, const double*,  std::vector<int>&, std::vector<unsigned int>&);
#endif



    //! Calculate the pair energy between two particles.
    /*! \param particle1
            The index of the first particle.

        \param position1
            The position vector of the first particle.

        \param orientation1
            The orientation vector of the first particle.

        \param status1
            The status vector of the patches in particle 1

        \param particle2
            The index of the second particle.

        \param position2
            The position vector of the second particle.

        \param orientation2
            The orientation vector of the second particle.

        \param status2
            The status vector of the patches in particle 2
        \return
            The pair energy between particles 1 and 2.
     */
#ifndef ISOTROPIC
    virtual double computePairEnergy(unsigned int, const double*, const double*, const unsigned int*, unsigned int, const double*, const double*, const unsigned int*);
#else
    virtual double computePairEnergy(unsigned int, const double*, unsigned int, const double*);
#endif

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
#ifndef ISOTROPIC
    virtual unsigned int computeInteractions(unsigned int, const double*, const double*, unsigned int*);
#else
    virtual unsigned int computeInteractions(unsigned int, const double*, unsigned int*);
#endif

    //! Apply any post-move updates for a given particle.
    /*! \param particle
            The particle index.

        \param position
            The position of the particle following the virtual move.

        \param orientation
            The orientation of the particle following the virtual move.
        
        \param status
            The status of the patches particle following the virtual move
    */
#ifndef ISOTROPIC
    virtual void applyPostMoveUpdates(unsigned int, const double*, const double*, const unsigned int*);
#else
    virtual void applyPostMoveUpdates(unsigned int, const double*);
#endif

    //! Get the average pair energy.
    /*! \return
            The average pair energy.
     */
    virtual double getEnergy();
    unsigned int maxInteractions;	//!< The maximum number of interactions per particle.

protected:
    double interactionEnergy;           //!< Interaction energy scale (in units of kBT).
    double interactionRange;            //!< Size of interaction range (in units of particle diameter).
    double squaredCutOffDistance;       //!< The squared cut-off distance.
};

#endif  /* _MODEL_H */
