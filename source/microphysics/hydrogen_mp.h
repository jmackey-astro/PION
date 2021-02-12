///
/// \file hydrogen_mp.h
/// \author Jonathan Mackey
/// \date 08.03.2011
///
/// This contains collected classes and functions for hydrogen microphysics.
///
/// Modifications:
/// - 2011.03.08 JM: written.
/// - 2011.03.14 JM: Added photoionisation inheritance.  removed redundant PI
/// x-section.
/// - 2011.03.29 JM: Added Hi_coll_ion_rates() function.

#ifndef HYDROGEN_MP_H
#define HYDROGEN_MP_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "hydrogen_photoion.h"
#include "hydrogen_recomb_Hummer94.h"

///
/// Hydrogen chemistry class, bringing together a lot of functions which
/// calculate heating/cooling rates and ionisation/recombination rates for
/// atomic Hydrogen.
///
class Hydrogen_chem : public Hummer94_Hrecomb, public hydrogen_photoion {
  public:
    Hydrogen_chem();
    ~Hydrogen_chem();

    //
    // Public functions inherited from Hummer94_Hrecomb:
    //
    // Hii_rad_recomb_rate(T)
    // Hii_rad_recomb_cooling(T)
    // Hii_total_cooling(T) (case B recomb +Bremsstrahlung)
    //
    // Public functions inherited from hydrogen_photoion
    // void Setup_photoionisation_rate_table(Tstar,Rstar,etc) (cgs)
    // double Hi_multifreq_photoionisation_rate(NH,nH,Vshell) (cgs)
    // double Hi_multifreq_photoionisation_heating_rate(NH,nH,Vshell) (cgs)
    // double Hi_monochromatic_photo_ion_xsection(E) (cgs)
    //

    ///
    /// Collisional ionisation rate, using the fitting formulae from
    /// Voronov (1997) ADANDT, 65, 1. (Atomic Data And Nuclear Data Tables)
    /// The rate is returned in cm^3/s, and should be multiplied by n_{e}*n_{HI}
    /// to get a rate per unit volume.
    ///
    double Hi_coll_ion_rate(double  ///< electron temperature
    );

    ///
    /// Collisional ionisation cooling rate (calculated as for the ionisation
    /// rate, but assuming that the binding energy is lost from the system in
    /// every ionisation, assumed to be from the ground state).
    ///
    double Hi_coll_ion_cooling_rate(double  ///< electron temperature
    );

    ///
    /// Returns the collisional ionisation rate and cooling rate together
    /// (saves computation).
    ///
    void Hi_coll_ion_rates(
        double,   ///< electron temperature
        double*,  ///< ionisation rate
        double*   ///< cooling rate
    );

    ///
    /// Collisional excitation cooling rate, for a given electron temperature.
    /// Result is returned in units of erg.cm^3/s and should be multiplied by
    /// n_{e}*n_{HI} to give a cooling rate per unit volume.
    /// Tabulated values from Raga, Mellema, \& Lundqvist, (1997,ApJS,109,517),
    /// which are based on data from Aggarwal (1993).
    ///
    double Hi_coll_excitation_cooling_rate(double  ///< electron temperature
    );

    ///
    /// Energy added to the gas per photoionisation.  This is just the energy
    /// of the photon minus the binding energy.  Returns the heating rate per
    /// photoionisation in ergs.
    ///
    double Hi_monochromatic_photo_ion_heating(const double  ///< photon energy
    );

  private:
    ///
    /// This function will set up the spline interpolation table, and slopes for
    /// extrapolation outside the range of validity of the table.
    ///
    void setup_Hi_coll_excitation_rate();

    int cx_Nspl;       ///< n_elements in spline fit for coll. ex.
    int cx_spline_id;  ///< id of spline in interpolation class
    double *cx_T,      ///< list of log10(T) values for collisional excitation.
        *cx_rate;  ///< list of log10(Rate) values for collisional excitation.
    double cx_MinSlope,  ///< logarithmic slope for extrapolation of collisional
                         ///< excitation table.
        cx_MaxSlope,     ///< logarithmic slope for extrapolation of collisional
                         ///< excitation table.
        cx_minT,  ///< minimum temperature for collisional excitation table.
        cx_maxT;  ///< maximum temperature for collisional excitation table.
};

#endif  // HYDROGEN_MP_H
