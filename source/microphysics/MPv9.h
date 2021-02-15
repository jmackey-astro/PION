///
/// \file MPv10.h
/// \author Maggie Goulden
/// \date 2018.10
///
/// modifications:
/// - 2018.10 .

#ifndef MPv10_H
#define MPv10_H

/// Description:
/// A multi-ion chemistry solver with non-equilibrium ionization and
/// recombination of different species, and with elemental abundances
/// that are input for each solve (i.e. can vary within simulation).
///
/// The integration method uses the CVODES solver from the SUNDIALS
/// package by (Cohen, S. D., & Hindmarsh, A. C. 1996, Computers in
/// Physics, 10, 138) available from
///  https://computation.llnl.gov/casc/sundials/main.html
/// The method is backwards differencing (i.e. implicit) with Newton
/// iteration.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "cooling_SD93_cie.h"
#include "cvode_integrator.h"
#include "hydrogen_mp.h"
#include "hydrogen_photoion.h"
#include "hydrogen_recomb_Hummer94.h"
#include "microphysics_base.h"
#include <cstring>
#include <iostream>
#include <set>
#include <sstream>
#include <vector>

/// this is the max change in x or E which is integrated with Euler integration.
/// A larger change will be integrated with backward-differencing and Newton
/// iteration, which is more accurate, stable, but expensive.
#define EULER_CUTOFF 0.05

#define JM_RELTOL                                                              \
    1.0e-4  ///< relative-error tolerance (actual error can be larger).
#define JM_MINNEU 1.0e-20  ///< minimum neutral fraction i care about.
#define JM_MINERG 1.0e-17  ///< Minimum internal energy density I care about.

///
/// Integrator for the microphysics equations for the non-equilibrium ion
/// fraction of various ions, and for the internal energy density.
///
class MPv10 :
    public Hydrogen_chem,
    public cooling_function_SD93CIE,
    public microphysics_base,
    public cvode_solver {
  public:
    ///
    /// Constructor
    ///
    MPv10(
        const int,              ///< grid dimensions
        const int,              ///< Coordinate System flag
        const int,              ///< Total number of variables in state vector
        const int,              ///< Number of tracer variables in state vector.
        const std::string*,     ///< List of what the tracer variables mean.
        struct which_physics*,  ///< pointer to extra physics flags.
        struct rad_sources*,    ///< radiation sources.
        const double            ///< EOS Gamma
    );

    /// MPv10::MPv10
    /// Destructor
    ///
    ~MPv10();

    ///
    /// Function for updating vectors according to species found in tracer list.
    ///
    void species_tracer_initialise(
        const std::string*,  ///< List of what the tracer variables mean.
        int,  ///< index of current tracer in for loop (i in s=tracers[i])
        std::string,  /// < current tracer in for loop
        std::string,  ///< element symbol, e.g. "He", "H"
        int,  ///< element length, e.g. "H" is of length 1, "He" of length 2.
        int,  ///< element index in N_species_by_elem, used in for loops for
              ///< densities etc
        int   ///< length of tracers vector
    );

    ///
    /// The NON-RT MICROPHYSICS update function.
    ///
    /// A primitive vector is input, and lots of optional extra data in
    /// the last argument, and the ion fraction(s) and internal energy
    /// are updated by the requested timestep.  Results are saved to
    /// the destination vector, which can be the same pointer as the
    /// input vector.
    ///
    int TimeUpdateMP(
        const pion_flt*,  ///< Primitive Vector to be updated.
        pion_flt*,        ///< Destination Vector for updated values.
        const double,     ///< Time Step to advance by.
        const double,     ///< EOS gamma.
        const int,        ///< Switch for what type of integration to use.
        double*           ///< Vector of extra data (column densities, etc.).
    );

    ///
    /// UNUSED FUNCTION!!
    ///
    int TimeUpdate_RTsinglesrc(
        const pion_flt*,  ///< Primitive Vector to be updated.
        pion_flt*,        ///< Destination Vector for updated values.
        const double,     ///< Time Step to advance by.
        const double,     ///< EOS gamma.
        const int,        ///< Switch for what type of integration to use.
        const double,     ///< flux in per unit length along ray (F/ds or L/dV)
        const double,     ///< path length ds through cell.
        const double,     ///< Optical depth to entry point of ray into cell.
        double*  ///< return optical depth through cell in this variable.
    )
    {
        cout << "MPv10::TimeUpdate_RTsinglesrc() not implemented\n";
        return 1;
    }

    ///
    /// Not used for this class, so far
    ///
    virtual int TimeUpdateMP_RTnew(
        const pion_flt*,  ///< Primitive Vector to be updated.
        const int,        ///< Number of UV heating sources.
        const std::vector<struct rt_source_data>&,
        ///< list of UV-heating column densities and source properties.
        const int,  ///< number of ionising radiation sources.
        const std::vector<struct rt_source_data>&,
        ///< list of ionising src column densities and source properties.
        pion_flt*,     ///< Destination Vector for updated values
        const double,  ///< Time Step to advance by.
        const double,  ///< EOS gamma.
        const int,     ///< Switch for what type of integration to use.
        double*        ///< any returned data (final temperature?).
    )
    {
        cout << "MPv10::TimeUpdateMP_RTnew() not implemented\n";
        return 1;
    }

    ///
    /// Returns the gas temperature.  This is only needed for data output, so
    /// there is no need to make it highly optimized.
    /// - Is threadsafe.
    ///
    double Temperature(
        const pion_flt*,  ///< primitive vector
        const double      ///< eos gamma
    );

    ///
    /// Set the gas temperature to a specified value.
    /// Needed if you want this feature to set initial conditions.
    ///
    int Set_Temp(
        pion_flt*,     ///< primitive vector.
        const double,  ///< temperature
        const double   ///< eos gamma.
    );

    ///
    /// This returns the minimum timescale of the times flagged in the
    /// arguments.
    ///
    virtual double timescales(
        const pion_flt*,  ///< Current cell.
        const double,     ///< EOS gamma.
        const bool,       ///< set to 'true' if including cooling time.
        const bool,       ///< set to 'true' if including recombination time.
        const bool        ///< set to 'true' if including photo-ionsation time.
    );

    ///
    /// This returns the minimum timescale of all microphysical
    /// processes, including reaction times for each species and the
    /// total heating/cooling time for the gas.
    /// It requires the radiation field as an input, but this is not
    /// used for now so the vectors can be null.
    ///
    virtual double timescales_RT(
        const pion_flt*,  ///< Current cell.
        const int,        ///< Number of UV heating sources.
        const std::vector<struct rt_source_data>&,
        ///< list of UV-heating column densities and source properties.
        const int,  ///< number of ionising radiation sources.
        const std::vector<struct rt_source_data>&,
        ///< list of ionising src column densities and source properties.
        const double  ///< EOS gamma.
    );

    ///
    /// Initialise microphysics ionisation fractions to an equilibrium
    /// value. This is optionally used in the initial condition
    /// generator.  Not implemented here.
    ///
    int Init_ionfractions(
        pion_flt*,     ///< Primitive vector to be updated.
        const double,  ///< eos gamma.
        const double   ///< optional gas temperature to end up at.
    )
    {
        cout << "MPv10::Init_ionfractions() not implemented\n";
        return 1;
    }

    ///
    /// Return index of tracer for a given string.
    ///
    int Tr(const string  ///< name of tracer we are looking for.
    );

    ///
    /// Get the total cooling rate.  This is for postprocessing the
    /// simulation data only -- IT IS NOT OPTIMISED FOR SPEED.
    ///
    virtual double total_cooling_rate(
        const pion_flt*,  ///< Current cell values.
        const int,        ///< Number of UV heating sources.
        const std::vector<struct rt_source_data>&,
        ///< list of UV-heating column densities and source properties.
        const int,  ///< number of ionising radiation sources.
        const std::vector<struct rt_source_data>&,
        ///< list of ionising src column densities and source properties.
        const double  ///< EOS gamma.
    )
    {
        cout << "MPv10::total_cooling_rate() not implemented\n";
        return 1;
    }

    ///
    /// Get the total recombination rate for an ion, given the input
    /// state vector.
    ///
    virtual double get_recombination_rate(
        const int,        ///< ion index in tracer array (optional).
        const pion_flt*,  ///< input state vector (primitive).
        const double      ///< EOS gamma (optional)
    )
    {
        cout << "MPv10::get_recombination_rate() not implemented\n";
        return 1;
    }

    ///
    /// Return the H mass fraction
    ///
    virtual inline double get_X_H() { return EP->H_MassFrac; }

  protected:
    ///
    /// convert state vector from grid cell into local microphysics vector.
    ///
    virtual int convert_prim2local(
        const pion_flt*,  ///< primitive vector from grid cell (length nv_prim)
        double*           ///< local vector [x(H0),E](n+1).
    );

    ///
    /// Convert local microphysics vector into state vector for grid cell.
    /// This is the inverse of convert_prim2local.
    ///
    virtual int convert_local2prim(
        const double*,    ///< local (updated) vector [x(H0),E](n+1).
        const pion_flt*,  ///< input primitive vector from grid cell (length
                          ///< nv_prim)
        pion_flt*  ///< updated primitive vector for grid cell (length nv_prim)
    );

    ///
    /// returns gas temperature according to E=nkT/(g-1) with n=1.1*nH*(1+x_in),
    /// appropriate for a gas with 10% Helium by number, and if He is singly
    /// ionised whenever H is.
    ///
    virtual double get_temperature(
        double*,  // const double, ///< nH //< y_ion_fraction (by y_ion_index)
        vector<pion_flt>&,  ///< y_ion_number_density
        const double        ///< E_int (per unit volume)
    );

    ///
    /// Returns total number of particles.
    ///
    virtual double get_ntot(
        double*,  // const double, ///< nH //< y_ion_fraction (by y_ion_index)
        vector<pion_flt>&  ///< y_ion_number_density
    );

    ///
    /// Set the size of local vectors, and index them with integers.
    ///
    virtual void setup_local_vectors();

    //
    // ********* STUFF FROM THE mp_v2_aifa CLASS **********
    //
    N_Vector y_in,    ///< current y-vector
        y_out;        ///< output y-vector
    int N_equations;  ///< number of equations i.e. number of species in Y.
    int N_extradata;  ///< number of elements in user-data array.

    //
    // Any data which can be calculated/set at the start of the simulation
    // can be defined here.
    //
    double k_B;              ///<  Boltzmanns constant.
    double m_p;              ///< Mass of proton.
    double m_H;              ///< Mass of hydrogen (grams).
    double m_He;             ///< Mass of helium.
    const int ndim;          ///< Number of dimensions in grid.
    const int nv_prim;       ///< Number of variables in state vector.
    const double eos_gamma;  ///< EOS gamma for ideal gas.
    const int coord_sys;     ///< Coordinate System flag
    double gamma_minus_one;  ///< as named.
    double Min_NeutralFrac;  ///< minimum H0 fraction allowed (eps=1.0e-12)
    double Max_NeutralFrac;  ///< Maximum H0 fraction allowed (1-eps)
    double mean_mass_per_H;  ///< mean mass per hydrogen nucleon, should be
                             ///< about 2.34e-24;
    double JM_NELEC;         ///< Number of electrons per ionised H atom.
    double JM_NION;          ///< Number of ions per ionised H atom.
    double METALLICITY;      ///< Metallicity of gas, in units of solar.

    int nvl;                    ///< number of variables in local state vector.
    int lv_eint;                ///< internal energy local variable index.
    int lv_H0;                  ///< neutral hydrogeen local variable index.
    int lv_y_ion_index_offset;  ///< gives the index at which ions first occur
                                ///< in primitive vector, maps to first index of
                                ///< local vector

    int N_elem;
    int N_species;
    int N_eqns;       ///< := n species
    int N_cons_eqns;  /// < := n_elem + 1, due to E_int.

    /// ===========================================================================
    ///               Vectors to Access Primitive / Local vectors
    /// ===========================================================================
    std::vector<int> X_mass_frac_index;  /// < primitive vector indices, used to
                                         /// trace X_H etc, like pv_Hp.
    std::vector<int>
        y_ion_index_prim;  ///< index matching y_ion mass fraction in
                           ///< prim vector, analogous to pv_Hp before.
    std::vector<int>
        y_ion_index_local;  ///< index matching y_ion fraction in local vector.

    std::vector<int> H_ion_index;   ///< Locates position of ion with N+1
                                    ///< electrons missing, e.g. H_ion_index[0]
                                    ///< -> H+ position. Used with MPv10::Tr().
    std::vector<int> He_ion_index;  ///""

    /// ===========================================================================
    ///               Vectors to Access Adjacent Ions
    /// ===========================================================================
    std::vector<int>
        y_ip1_index_local;  ///< local index of ion p1 (plus 1); ion one stage
                            ///< higher. If ion doesn't exist, sets to -1.
    std::vector<int>
        y_im1_index_local;  ///< local index of ion m1 (minus 1); ion one stage
                            ///< lower. if lower stage ion is neutral species,
                            ///< sets to -2.
                            ///   also, if we're just not tracking lower stage
                            ///   (i.e. "doesn't exist"), also sets to -1.
    std::vector<int>
        y_ip1_index_tables;  ///< index of the ip1 species in tables. This
                             ///< follows the order "H0", "H1+", "He0", "He1+",
                             ///< "He2+" etc.
    std::vector<int> y_ion_index_tables;
    std::vector<int>
        y_im1_index_tables;  ///< index of the im1 species in tables.
                             ///< Follows same order as above.

    /// ===========================================================================
    ///           Vectors to Store Other Ion / Element info
    /// ===========================================================================
    std::vector<int>
        y_ion_num_elec;  ///< index matching number of electrons corresponding
                         ///< to each y_ion (e.g. C6+ = 6)
    std::vector<pion_flt>
        X_elem_atomic_mass;  /// < vector of atomic masses corresponding to
                             /// X_elem_index, e.g. for ["H", "He"],
                             /// X_atom_mass=[1.6738e-24, 6.6464764e-24
    std::vector<pion_flt> X_elem_number_density;
    std::vector<int> N_species_by_elem;  ///< records # species in each element,
                                         ///< to iterate over later.

    /// ===========================================================================
    ///  Variables to store table (temp, ionisation / recomb rates,
    ///  etc)information
    /// ===========================================================================
    std::vector<pion_flt> Temp_Table;
    std::vector<float> recomb_rate_table;
    std::vector<pion_flt> ionise_rate_table;
    std::vector<pion_flt> recomb_slope_table;
    std::vector<pion_flt> ionise_slope_table;

    std::vector<pion_flt>
        ionisation_potentials;  // array of ionisation potentials, in the same
                                // order as ionise_slope

    const float T_min;
    const float T_max;
    const int
        Num_temps;  // NB this needs to be const so can initialise arrays
                    // with it. If this is >=1e4, get stack overflow errors.
    float delta_log_temp;

    int pv_H1p;  ///< legacy, should ideally remove

    int N_diff_srcs,   ///< No diffuse sources --> 0, otherwise --> 1
        N_ion_srcs,    ///< No ionising sources --> 0, otherwise --> 1
        ion_src_type;  ///< Either RT_EFFECT_MFION or RT_EFFECT_PION_MONO.

    //---------------------------------------------------------------------------
    //-------------- FUNCTIONS DERIVED FROM BASE CLASS FOLLOW
    //-------------------
    //---------------------------------------------------------------------------
  public:
    ///
    /// calculate dy/dt for the vector of y-values.
    ///
    virtual int ydot(
        double,  ///< current time (probably not needed for rate equations)
        const N_Vector,  ///< current Y-value
        N_Vector,        ///< vector for Y-dot values
        const double*    ///< extra user-data vector, P, for evaluating
                         ///< ydot(y,t,p)
    );

    ///
    /// Calculate the Jacobian matrix d(dy_i/dt)/dy_j for a vector of y-values.
    ///
    virtual int Jacobian(
#if defined CVODE2
        int,  ///< N (not sure what this is for!)
#endif
        double,          ///< time, t
        const N_Vector,  ///< current Y-value
        const N_Vector,  ///< vector for Y-dot values
        const double*,   ///< extra user-data vector, P, for evaluating
                         ///< ydot(y,t,p)
        CVMatrix         ///< Jacobian matrix
    )
    {
        cout << "Jacobian not implemented in MPv10!\n";
        return 1;
    }

    ///
    /// Get the number of extra parameters and the number of equations.
    ///
    virtual void get_problem_size(
        int*,  ///< number of equations
        int*   ///< number of parameters in user_data vector.
    );

  protected:
    ///
    /// set the relative and absolute error tolerances
    ///
    virtual void get_error_tolerances(
        double*,  ///< relative error tolerance (single value)
        double[]  ///< absolute error tolerance (array)
    );

    //---------------------------------------------------------------------------
    //-------------- END OF FUNCTIONS DERIVED FROM BASE CLASS
    //-------------------
    //---------------------------------------------------------------------------

    //
    // Constant data in the cell, received from cell data.
    //
    double mpv_nH,   ///< total hydrogen number density at current cell.
        mpv_Vshell,  ///< geometric factor in point source flux calculation
                     ///< (\sim 4\pi R^2 dR).
        mpv_Tau0,   ///< Optical depth of neutral hydrogen to front edge of cell
                    ///< at 13.6eV
        mpv_dTau0,  ///< Optical depth of neutral hydrogen through cell
                    ///< (=nH*dS*(1-x)*sigma) at 13.6eV.
        mpv_G0_UV,  ///< UV heating flux, including attenuation, F*exp(-1.9Av).
        mpv_G0_IR,  ///< Heating due to UV flux re-radiated in IR and
                    ///< re-absorbed, F*exp(-0.05Av).
        mpv_NIdot,  ///< photon luminosity of monochromatic ionising source
                    ///< (ionising photons/s).
        mpv_delta_S;  ///< path length through current cell.
};

#endif  // MPv10_H
