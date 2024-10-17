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

#ifndef _SINGLEPARTICLEMOVE_H
#define _SINGLEPARTICLEMOVE_H

#include <vector>

#include "MersenneTwister.h"


/*! \file SingleParticleMove.h
    \brief A class for executing single particle Move Monte Carlo moves (translations and rotations).
*/

// FORWARD DECLARATIONS

class  Model;
struct Particle;
class  MersenneTwister;

//! Container for storing move parameters.
struct MoveParams
{
    unsigned int seed;                          //!< Index of the seed particle.
    bool isRotation;                            //!< Whether the move is a rotation.
    bool isChange;                              //!< Whether the move is a change of state.
    double stepSize;                            //!< The magnitude of the trial move.
    std::vector<double> trialVector;            //!< Vector for trial move.
    Particle preMoveParticle;                   //!< Particle state before the trial move.
    Particle postMoveParticle;                   //!< Particle state before the trial move.

    std::vector<unsigned int> preMoveBonds; //!< the number of particles the neighbours are in contact with   
    std::vector<unsigned int> postMoveBonds; //!< the number of particles the neighbours are in contact with   

};

class SingleParticleMove
{
public:
    //! Constructor.
    /*! \param model_
            A pointer to the model object.
        \param maxTrialTranslation_
            The maximum trial translation (in units of the reference particle diameter).
        \param maxTrialRotation_
            The maximum trial rotation.
        \param probTranslate_
            The probability of performing a translation move (versus a rotation).
        \parm probChange
            The probability of changing the status of a patch (versus rotation and translation
        \param isIsotropic_
            Whether the potential is isotropic.
     */
    SingleParticleMove(MersenneTwister&, Model*, double, double, double, double, double, bool);

    //! Overloaded ++ operator. Perform a single step.
    void operator ++ (const int);

    //! Overloaded += operator. Perform "n" steps.
    void operator += (const int);

    //! Perform a single trial move.
    void step();

    //! Perform a specified number of trial moves.
    /*! \param nSteps
            The number of attempted trial moves.
     */
    void step(const int);

    //! Get the number of attempted moves.
    /*! \return
            The number of attempted moves.
     */
    unsigned long long getAttempts() const;

    //! Get the number of accepted moves.
    /*! \return
            The number of accepted virtual moves.
     */
    unsigned long long getAccepts() const;

    //! Get the number of accepted rotation moves.
    /*! \return
            The number of accepted rotation moves.
     */
    unsigned long long getRotations() const;

     //! Get the number of accepted changes move
    /*! \return
            The number of accepted changes moves.
    */
    unsigned long long getChanges() const;

    //! Reset statistics.
    void reset();

    MersenneTwister& rng;                        //!< Random number generator.

private:
    MoveParams moveParams;                      //!< Parameters for the trial move.
    Model* model;                               //!< A pointer to the model object
    double surf_interaction;                    //!< The interaction with the surface
    unsigned long long nAttempts;               //!< Number of attempted moves.
    unsigned long long nAccepts;                //!< Number of accepted moves.
    unsigned long long nRotations;              //!< Number of accepted rotations.
    unsigned long long nChanges;                //!< Number of accepted changes
    double maxTrialTranslation;                 //!< The maximum trial translation (in units of the reference diameter).
    double maxTrialRotation;                    //!< The maximum trial rotation.
    double probTranslate;                       //!< The relative probability of translational moves (vs rotations).
    double probChange;                          //!< The relative probability of Changes vs (movement)
    bool is3D;                                  //!< Whether the simulation is three-dimensional.
    bool isIsotropic;                           //!< Whether the potential is isotropic.
    double energyChange;                        //!< Energy change resulting from trial move.
    double initialEnergy;                       //!< Initial interaction energy of state a
    //! Propose a trial particle translation/rotation.
    void proposeMove();
    void proposeMove_post();

    //! Determine whether move is accepted.
    bool accept();

    //! Calculate an unbiased rotation vector in 3D (Beard & Schlick, BJ 85 2973 (2003)).
    /*! \param v1
            The vector about which to rotate (either the position or orientation).

        \param v2
            Rotation unit vector.

        \param v3
            The rotation vector.

        \param angle
            Trial rotation angle.
     */
    void rotate3D(std::vector<double>&, std::vector<double>&, std::vector<double>&, double);

    //! Calculate a simple in plane rotatation vector.
    /*! \param v1
            The vector to rotate (either the position or orientation).

        \param v2
            The rotation vector.

        \param angle
            Trial rotation angle.
     */
    void rotate2D(std::vector<double>&, std::vector<double>&, double);

    //compute a possible new state for a patch and compute free energy of change
    /*! \param v1
            The id of the patches to change
        \param v2
            The proability register for each state
        \param state_probability to switch to vector in case of multiple possible states to change ad same time
            The probability to obtain a close state
    */
    double changepatch(std::vector<unsigned int>&, std::vector<double>&, double);

    //! Compute the norm of a vector.
    /*! \param vec
            A vector to normalize
        \return
            The norm of the vector.
     */
    double computeNorm(std::vector<double>&);
};

#endif  /* _SINGLEPARTICLEMOVE_H */
