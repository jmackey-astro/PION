/// \file riemannMHD.h
/// \brief Class definition for linear MHD Riemann Solver.
///
///  \author Jonathan Mackey
///
/// Written 2006-11-20
/// Modified :
///  - 2007-07-13 added in new function for SetDirection() and SetAvgState()
///  - 2007-07-16 fixed problem in solver2codevars() introduced by
///  SetDirection() on 13/7.
///  - 2009-10-24 New class structure, so that it inherits from the MHD
///  equations class.
///
/// This is a linear MHD Riemann Solver for the ideal MHD equations.
/// It is a Roe-type solver, linearising the Jacobian matrix with an
/// average state and jumping across waves from the left and/or right
/// to get to the state at x=0.
///
/// - 2010.12.23 JM: updated so that there is no riemann_base class,
///   and so that RS_left/right/pstar/meanp are all private vars.
///   Removed SetAvgState() (moved to eqns_mhd_ideal).
/// - 2015.06.07 JM: Tidied up a bit!
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)

#ifndef RIEMANNMHD_H
#define RIEMANNMHD_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "equations/eqns_mhd_adiabatic.h"
#include "riemann.h"

#include <string>

/// \brief MHD Waves Enum */
enum waves {
    FN = 0,  ///< Negative Fast Magnetic Wave
    AN = 1,  ///< Negative Alfven Wave
    SN = 2,  ///< Negative Slow Magnetic Wave
    CT = 3,  ///< Contact Discontinuity
    SP = 4,  ///< Positive Slow Magnetic Wave
    AP = 5,  ///< Positive Alfven Wave
    FP = 6   ///< Positive Fast Magnetic Wave
};

/// \brief MHD Riemann Solver Primitive Variables Vector
///
/// This vector has no B_x, because that must be kept constant in the
/// 1D Riemann Problem. (The fact that it is not constant leads to
/// errors and monopoles, but that's to be sorted out elsewhwere).
/// So this enum is to be used only in the Riemann Solver where I have
/// seven variables instead of eight.
///
/// Using this, I keep the B_x info in the eigth element, but for things
/// like dot products I can just use the first seven elements and put it
/// in a for loop.
///
enum rsvars {
    RRO = 0,  ///< Density
    RPG = 1,  ///< Gas pressure
    RVX = 2,  ///< x-velocity
    RVY = 3,  ///< y-velocity
    RVZ = 4,  ///< z-velocity
    RBY = 5,  ///< y-magnetic field
    RBZ = 6,  ///< z-magnetic field
    RBX = 7   ///< x-magnetic field (not needed except as a parameter)
};

/// \brief Adiabatic MHD linear Riemann Solver.
///
/// \author Jonathan Mackey
///
/// In this class, the boundary between left and right states is at x=0, and
/// is in the y-z plane.  The left state is x<0, and the right state is x>0.
///
/// \section References
/// Falle, Komissarov & Joarder, 1998, MNRAS, 297, 265.\n
/// Roe \& Balsara, 1996, SIAM J. Applied Math., 56, 57.
///
///
class riemann_MHD :
    // virtual public riemann_base,
    virtual public eqns_mhd_ideal {
  public:
    /// \brief Constructor
    ///
    /// This sets up memory for the various state vectors needed for the solver.
    /// Also receives mean values of problems it is likely to encounter, for
    /// reference purposes, when deciding whether a quantity is 'small' or not.
    ///
    riemann_MHD(
        const int,        ///< Length of state vectors, nvar.
        const pion_flt*,  ///< Mean values of primitive variables on grid
                          ///< [vector, length nvar]
        const double      ///< Gamma for state vector.
    );

    /// \brief Destructor: deletes dynamically allocated arrays.
    ~riemann_MHD();

    /// Public solve function: called to solve the riemann problem
    ///
    /// This function is the only public interface to the solver.  It takes in
    /// pointers to the left and right state vectors, and a pointer to the
    /// resulting state vector.
    ///
    /// \retval 0 success
    /// \retval 1 fatal failure in method, so should bug out and exit code.
    /// \retval 2 non-fatal failure, e.g. got negative pressure so return vacuum
    /// conditions.
    ///
    int JMs_riemann_solve(
        const pion_flt*,  ///< Left state vector.
        const pion_flt*,  ///< Right state vector.
        pion_flt*,        ///< Result state vector.
        const int,        ///< mode to solve (FLUX_RSlinear is only option!)
        const double      ///< Gas eos constant gamma.
    );

  private:
    int RS_mode;  ///< Mode of solution, only mode=1 works so far, for the
                  ///< linear solver.

    int RS_nvar;  ///< The number of active variables/waves in the Riemann
                  ///< problem. Hard-coded to seven for simplicity.

    double ansBX;  ///< B_x is just a parameter, so set it at the start of the
                   ///< sovle function, and assign to this variable.

    double ch,   ///< Hydrodynamical sound speed (in average state).
        ca,      ///< Alfven speed (in average state).
        cs,      ///< Slow magnetic wave speed (in average state).
        cf,      ///< Fast magnetic wave speed (in average state).
        bx,      ///< Normalised x-component of magnetic
                 ///< field,\f$=B_x/\sqrt{\rho}\f$.
        bt,      ///< Normalised tangential component of magnetic
                 ///< field,\f$=\sqrt{(B_y^2+B_z^2)/\rho}\f$.
        betay,   ///< This is \f$B_y/\sqrt{B_y^2+B_z^2}\f$.
        betaz,   ///< This is \f$B_z/\sqrt{B_y^2+B_z^2}\f$.
        alphaf,  ///< The fast speed normalisation constant
                 ///< \f$=(a^2-c_s^2)/(c_f^2-c_s^2)\f$.
        alphas;  ///< The slow speed normalisation constant
                 ///< \f$=(c_f^2-a^2)/(c_f^2-c_s^2)\f$.

    pion_flt* RS_pdiff;     ///< Difference between left and right states (for
                            ///< calculation).
    pion_flt* RS_evalue;    ///< The eignevalues of the matrix \f$\bar{A}\f$.
    pion_flt* RS_strength;  ///< The wavestrengths for each wave.

    /// Left state vector (local copy)
    pion_flt* RS_left;
    /// Right state vector (local copy)
    pion_flt* RS_right;
    /// Resolved state vector (local copy)
    pion_flt* RS_pstar;
    /// Mean state vector
    pion_flt* RS_meanp;

    /// The elements of the 7 left eigenvectors... evec[E-val][Element]
    pion_flt RS_leftevec[7][7];
    /// The elements of the 7 right eigenvectors...evec[E-val][Element]
    pion_flt RS_rightevec[7][7];

    long int
        onaxis,     ///< Debugging counters for checking eigenvectors are right.
        offaxis,    ///< Debugging counters for checking eigenvectors are right.
        samestate,  ///< Debugging counter, for counting how many solves had the
                    ///< same left and right state.
        totalsolve;  ///< Counter for total number of Riemann Problems solved.

    ///
    /// Small number slightly larger than roundoff error for double
    /// precision. Used to make numbers positive when they should be,
    /// but aren't b/c of roundoff. For my machine, double precision
    /// accuracy is about 2e-16, so set this to 1.e-15.
    ///
    double smallB;

    /// Very small number, used to check that B_t is not zero.
    double tinyB;

    /// \brief assigns the data passed to the solver, to class member variables.
    ///
    /// Make sure the state vectors are in the format expected by the riemann
    /// solver!!!
    ///
    void assign_data(
        const pion_flt*,  ///< Pointer to left state vector
        const pion_flt*   ///< Pointer to right state vector
    );

    /// \brief  Calculates the average state vector, e.g. (P_L+P_R)/2
    ///
    /// This could have a number of ways of calculating an average state,
    /// but for now I just do a straight mean of the left and right states.
    ///
    void get_average_state();

    /// \brief Calculates the sound speeds for the left and right states.
    ///
    /// Calculates sound speeds (c_h, c_a, c_s, c_f)
    ///
    int get_sound_speeds();

    /// \brief  Calculates the eignevalues of the matrix \f$\bar{A}\f$
    ///
    /// These are the seven wavespeeds.  For a linear solver they are not the
    /// real wavespeeds, but speeds calculated from the average state which are
    /// assumed to apply everywhere.
    ///
    void get_eigenvalues();

    /// \brief  Check and output the values of the evectors and their dot
    /// products
    ///
    ///
    int check_evectors();

    /// \brief  This constructs evectors using some idea I came up with
    ///
    /// The vectors are from Falle et al., but I tried to be a bit smarter with
    /// the normalisation.  I can't remember, but I don't think it worked any
    /// better than Falle's norm.
    ///
    int my_evectors();

    /// \brief  This constructs evectors using the Roe and Balsara
    /// normalisation.
    ///
    /// Mostly Roe and Balsara norm, but the Alfven waves are still the same as
    /// in Falle et al.
    ///
    int RoeBalsara_evectors();

    /// \brief Evaluates dot product to two vectors
    ///
    /// Pass in the two vectors and their length.
    ///
    /// \retval value of dot product if successful.
    ///
    double dot_product(
        pion_flt*,  ///< Pointer to Vector 1.
        pion_flt*,  ///< Pointer to Vector 2.
        int         ///< Length of vectors.
    );

    /// \brief Calculates the wave strengths alpha_i
    ///
    /// The equation is easily derived, and is given in Falle et al., just after
    /// eq.8 as
    /// \f$ \alpha_i = \frac{\mathbf{l}^{[i]}\cdot\mathbf{(P_R-P_L)}}
    ///           {\mathbf{l}^{[i]}\cdot \mathbf{r}^{[i]}}\;. \f$
    ///
    /// \retval 0 success
    /// \retval 1 failure
    ///
    void calculate_wave_strengths();

    /// \brief  Calculates pdiff[]= (P_L-P_R) for use in getting P*
    /// \retval 0 success
    /// \retval 1 failure
    ///
    void getPdiff();

    /// \brief Gets the appropriate value of P* from left and right and makes
    /// sure they match.
    ///
    ///
    /// \retval 0 success
    /// \retval 1 failure
    ///
    int get_pstar();

    /// \brief Checks the error code, and if it's non-zero, prints it and a
    /// message
    ///
    void failerror(
        int,         ///< Error code.
        std::string  ///< Error message
    );

    /// \brief Changes the order of variables so they are back to the code order
    ///
    /// The riemann solver orders the variables according to the enum rsvars,
    /// while the main code orders them according to the enum primitive, so this
    /// function changes the order of the B-field elements from rsvars to
    /// primitive.
    ///
    ///
    void solver2codevars(pion_flt*  ///< Vector to convert.
    );

    /// \brief Change order of variables in state vector from code to solver
    /// variables.
    void code2solvervars(pion_flt*  ///< Vector to convert.
    );
};

#endif
