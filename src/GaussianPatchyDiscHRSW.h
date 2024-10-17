#ifndef _GAUSSIANPATCHYDISCHRSW_H
#define _GAUSSIANPATCHYDISCHRSW_H

#include "GaussianPatchyDisc.h"


/*! \file GaussianPatchyDiscHRSW.h
*/

//! Class defining the Patchy-Disc-HRSW potential.
class GaussianPatchyDiscHRSW : public GaussianPatchyDisc
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
    GaussianPatchyDiscHRSW(Box&, std::vector<Particle>&, CellList&, Top&, std::vector<unsigned int>&, unsigned int, double, double);
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
    virtual double computePairEnergy(unsigned int, const double*, const double*, unsigned int, const double*, const double*);

};

#endif  /* _GAUSSIANPATCHYDISCHR_H */
