///
/// \file mpv9_HHe.h
/// \author Jonathan Mackey
/// \date 2013.08.19
///
/// This class is for running photoionisation calculations where it 
/// important to track the ionisation of Helium and Hydrogen.
///
///
/// Modifications:
/// - getting it written: mods up until 2013.08.23


#ifndef MPV9_HHE_H
#define MPV9_HHE_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#ifndef EXCLUDE_MPV9

#include "microphysics/microphysics_base.h"
#include "HHe_photoion.h"
#include "microphysics/cvode_integrator.h"


class mpv9_HHe :
  public MicroPhysicsBase,
  public HHe_photoion,
  public cvode_solver  
{
  public:
  ///
  /// Constructor
  ///
  mpv9_HHe(
          const int,              ///< Total number of variables in state vector
	  const int,              ///< Number of tracer variables in state vector.
	  const std::string &,    ///< List of what the tracer variables mean.
          struct which_physics *, ///< extra physics stuff.
          const double            ///< EOS gamma
	  );

  ///
  /// Destructor
  ///
  ~mpv9_HHe();


  //---------------------------------------------------------------------------
  //-------------- FUNCTIONS DERIVED FROM BASE CLASS FOLLOW -------------------
  //---------------------------------------------------------------------------
  ///
  /// This takes a copy of the primitive vector and advances it in time over
  /// the step requested, and at the end copies the updated vector into the
  /// destination vector.  For fully local microphysics but WITH radiative transfer,
  /// where the column densities for diffuse and direct radiation are included as 
  /// parameters.  The input list of column densities is ordered in a list of UV
  /// heating sources and a list of ionising sources.
  ///
  int TimeUpdateMP_RTnew(
          const double *, ///< Primitive Vector to be updated.
          const int,      ///< Number of UV heating sources.
          const std::vector<struct rt_source_data> &,
          ///< list of UV-heating column densities and source properties.
          const int,      ///< number of ionising radiation sources.
          const std::vector<struct rt_source_data> &,
          ///< list of ionising src column densities and source properties.
          double *,       ///< Destination Vector for updated values
          ///< (can be same as first Vector.
          const double,   ///< Time Step to advance by.
          const double,   ///< EOS gamma.
          const int, ///< Switch for what type of integration to use.
          ///< (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
          double *    ///< final temperature (not strictly needed).
          );

  ///
  /// Get the total recombination rate for an ion, given the input
  /// state vector.
  ///
  double get_recombination_rate(
          const int,      ///< ion index in tracer array.
          const double *, ///< input state vector (primitive).
          const double    ///< EOS gamma
          );

  ///
  /// Return number density of a given element.
  ///
  virtual double get_n_el(
        const double *, ///< primitive state vector.
        const int       ///< integer identifier for the element.
        );

  ///
  /// Returns the gas temperature.  This is only needed for data output, so
  /// there is no need to make it highly optimized.
  ///
  double Temperature(
            const double *, ///< primitive vector
            const double    ///< eos gamma
            );

  ///
  /// Set the gas temperature to a specified value.
  /// Again only needed if you want this feature in the initial condition generator.
  ///
  int Set_Temp(
          double *,     ///< primitive vector.
          const double, ///< temperature
          const double  ///< eos gamma.
          );

  ///
  /// This returns the appropriate timescales for ionisation or
  /// recombination, and heating or cooling, depending on how the
  /// module is set up.
  ///
  virtual double timescales_RT(
        const double *, ///< Current cell.
        const int,      ///< Number of UV heating sources.
        const std::vector<struct rt_source_data> &,
        ///< list of UV-heating column densities and source properties.
        const int,      ///< number of ionising radiation sources.
        const std::vector<struct rt_source_data> &,
        ///< list of ionising src column densities and source properties.
        const double   ///< EOS gamma.
        );

  ///
  /// Return index of tracer for a given string. (only hydrogen!)
  ///
  int Tr(const string ///< name of tracer we are looking for.
	);

  ///
  /// Get the ionisation cross section for an atom/ion at its
  /// threshold frequency.
  ///
  double get_th_xsection(
        const int id  ///< integer identifier for the ion.
        ) {return xsection_th(id);};


  ///
  /// The NON-RT MICROPHYSICS update function.
  /// Should never be called with this class.
  ///
  int TimeUpdateMP(
        const double *, ///< Primitive Vector to be updated.
        double *,       ///< Destination Vector for updated values.
        const double,   ///< Time Step to advance by.
        const double,   ///< EOS gamma.
        const int,      ///< Switch for what type of integration to use.
        double *        ///< Vector of extra data (column densities, etc.).
        ) {return DONT_CALL_ME;}
  
  ///
  /// This is only for the implicit time-integrator.
  ///
  int TimeUpdate_RTsinglesrc(
        const double *, ///< Primitive Vector to be updated.
        double *,       ///< Destination Vector for updated values.
        const double,   ///< Time Step to advance by.
        const double,   ///< EOS gamma.
        const int,      ///< Switch for what type of integration to use.
        const double,   ///< flux in per unit length along ray (F/ds or L/dV)
        const double,   ///< path length ds through cell.
        const double,   ///< Optical depth to entry point of ray into cell.
        double *        ///< return optical depth through cell in this variable.
        ) {return DONT_CALL_ME;}

  ///
  /// This returns the minimum timescale of the times flagged in the
  /// arguments.  Not implemented here because it is only for sims
  /// without radiative transfer.
  ///
  virtual double timescales(
          const double *, ///< Current cell.
          const double,   ///< EOS gamma.
          const bool, ///< set to 'true' if including cooling time.
          const bool, ///< set to 'true' if including recombination time.
          const bool  ///< set to 'true' if including photo-ionsation time.
          ) {return DONT_CALL_ME;}


  ///
  /// Initialise microphysics ionisation fractions to an equilibrium
  /// value.  Not implemented here.
  ///
  int Init_ionfractions(
        double *, ///< Primitive vector to be updated.
        const double, ///< eos gamma.
        const double  ///< optional gas temperature to end up at.
        ) {return DONT_CALL_ME;}

  //---------------------------------------------------------------------------
  //-------------- END OF FUNCTIONS DERIVED FROM BASE CLASS -------------------
  //---------------------------------------------------------------------------

  //---------------------------------------------------------------------------
  //-------------- FUNCTIONS DERIVED FROM CVODE CLASS FOLLOW ------------------
  //---------------------------------------------------------------------------
  public:
  ///
  /// calculate dy/dt for the vector of y-values (NOT IMPLEMENTED HERE).
  ///
  int ydot(
          double,         ///< current time (probably not needed for rate equations)
          const N_Vector, ///< current Y-value
          N_Vector,       ///< vector for Y-dot values
          const double *  ///< extra user-data vector, P, for evaluating ydot(y,t,p)
          );

  ///
  /// Calculate the Jacobian matrix d(dy_i/dt)/dy_j for a vector of y-values.
  ///
  int Jacobian(
          int,            ///< N (not sure what this is for! Must be internal)
          double,         ///< time, t
          const N_Vector, ///< current Y-value
          const N_Vector, ///< vector for Y-dot values
          const double *, ///< extra user-data vector, P, for evaluating ydot(y,t,p)
          DlsMat          ///< Jacobian matrix
          ) {cout <<"Jacobian not implemented in mp_explicit_H!\n"; return 1;}

  ///
  /// Get the number of extra parameters and the number of equations.
  ///
  void get_problem_size(int *, ///< number of equations
                        int *  ///< number of parameters in user_data vector.
                        );

  protected:
  ///
  /// set the relative and absolute error tolerances
  ///
  void get_error_tolerances(
          double *, ///< relative error tolerance (single value)
          double *  ///< absolute error tolerance (array)
          );

  //---------------------------------------------------------------------------
  //-------------- END OF FUNCTIONS DERIVED FROM CVODE CLASS ------------------
  //---------------------------------------------------------------------------

  protected:

  ///
  /// convert primitive state vector into local microphysics vector.
  ///
  virtual int convert_prim2local(
            const double *, ///< primitive vector [nv_prim]
            double *        ///< local vector [nvl].
            );

  ///
  /// Convert local microphysics vector into primitive state vector.
  /// This is the inverse of convert_prim2local.
  ///
  virtual int convert_local2prim(
            const double *, ///< local (updated) vector [nvl].
            const double *, ///< input primitive vector [nv_prim]
            double *       ///< updated primitive vector [nv_prim]
            );

  ///
  /// Interpret the input radiation sources into local optical depths
  ///
  void interpret_radiation_data(
          const int,  ///< number of heating sources (UNUSED)
          const std::vector<struct rt_source_data> &, ///< heatingsrc
          const int,  ///< number of ionising sources (==1)
          const std::vector<struct rt_source_data> & ///< ionisingsrc
          );

  ///
  /// returns gas temperature according to ion fraction
  ///
  virtual double get_temperature(
    const double * ///< local vector [nvl]
    );

  ///
  /// Returns number density of particles.
  ///
  double get_ntot(
    const double * ///< local vector [nvl]
    );

  ///
  /// Returns number density of electrons.
  ///
  double get_ne(
    const double * ///< local vector [nvl]
    );

  N_Vector 
    y_in,  ///< current y-vector
    y_out; ///< output y-vector
  int N_equations; ///< number of equations i.e. number of species in Y.
  int N_extradata; ///< number of elements in user-data array.

  const double gamma;           ///< EOS gamma for ideal gas
  const double gamma_minus_one; ///< as named
  struct which_physics *EP; ///< flags for which physics we are using

  double Min_Nfrac; ///< minimum H0 fraction allowed
  double Max_Nfrac; ///< Maximum H0 fraction allowed

  double mean_mass_per_H; ///< mean mass per hydrogen nucleon
  double X_H;  ///< number fraction of H (=1 by definition)
  double X_HE; ///< Number fraction of He compared to H: n(He)/n(H)

  const size_t nv_prim; ///< Number of variables in state vector.
  size_t       nvl;     ///< number of variables in local vector.

  size_t lv_nH;   ///< number density of H nucleons, local var. index
  size_t lv_E; ///< internal energy local variable index.

  size_t lv_H0;   ///< neutral hydrogen local variable index.
  size_t lv_He0;   ///< neutral helium local variable index.
  size_t lv_He1;   ///< singly-ionised helium local variable index.

  size_t pv_H0;  ///< primitive vector index for H0 fraction.
  size_t pv_He0; ///< primitive vector index for He0 fraction.
  size_t pv_He1; ///< primitive vector index for He1 fraction.

  //
  // Local data for ydot (NOT THREADSAFE!)
  //
  double nH;  ///< number density of H nucleons.
  double *tau; ///< column densities: tau0, tau1, tau2, tauD.
  double *dtau; ///< column densities: dtau0, dtau1, dtau2.
  double Vshell; ///< Vshell in the finite difference photoion. rate.
  double dS;     ///< column density through cell.

  double metallicity; ///< gas metallicity in units of solar (0.0142)
};

#endif //  if not EXCLUDE_MPV9

#endif // MPV9_HHE_H



