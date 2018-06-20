/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DIFF_SGS_TKE
#define DIFF_SGS_TKE

#include "diff.h"

class Diff_sgs_tke : public Diff
{
    public:
        Diff_sgs_tke(Model*, Input*);
        ~Diff_sgs_tke();

        void exec();
        void exec_viscosity();

        void set_values(); ///< this is not empty anymore, used to initialise statistics ... SvdLinden
        void exec_stats(Mask*); ///< function for additional statistics

        unsigned long get_time_limit(unsigned long, double);
        double get_dn(double);

        double tPr;

        #ifdef USECUDA
        // GPU functions and variables
        void prepare_device();
        void clear_device();
        #endif

        // Empty functions, there are allowed to pass.
        //void set_values() {}

    private:
        template<bool>
        void calc_strain2(double*,
                          double*, double*, double*,
                          double*, double*,
                          double*, double*,
                          double*, double*, double*);

        // NOOT: deze is reeds aangepast: double* toegevoegd voor sgs_tke-veld
        template<bool>
        void calc_evisc(double*,
                        double*, double*, double*, double*, double*,
                        double*, double*, double*,
                        double*, double*,
                        double*, double*, double*,
                        double, double);

        template<bool>
        void calc_evisc_neutral(double*,
                                double*, double*, double*, double*,
                                double*, double*,
                                double*, double*,
                                double, double);

        template<bool>
        void diff_u(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
        template<bool>
        void diff_v(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

        void diff_w(double*, double*, double*, double*, double*, double*, double*, double*, double*);

        template<bool>
        void diff_c(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

        double calc_dnmul(double*, double*, double);

        void init_stats(); ///< initialize additional statistiscs
        void calc_prandtl(double*, double*, double*, double*);
        void set_minimum_sgs_tke(double*); ///< enforce minimum value of SGS TKE

        double cs;

        // Added functions specifically for prognostic equation for SGS TKE ///< SJA van der Linden, 7 May 2018
        template<bool>
        void diff_sgs_tke(double*, double*, double*, double*, double*, double*, double*, double*, double*);

        void calc_sgs_tke_shear_tend_2nd   (double*, double*, double*);
        void calc_sgs_tke_buoyancy_tend_2nd(double*, double*, double*, double*, double*);
        void calc_sgs_tke_dissipation_2nd  (double*, double*, double*, double*);

        // For now, 4th order is not implemented ///< SJA van der Linden, 7 May 2018
        void calc_sgs_tke_shear_tend_4th() {}
        void calc_sgs_tke_buoyancy_tend_4th(){}
        void calc_sgs_tke_dissipation_4th() {}

        // For now place additional constants for SGS-TKE model (Deardorff; 1973,1980) here ///< SJA van der Linden, 7 May 2018
        const double ap  = 1.5;
        const double cf  = 2.0;//2.5; temporary change as prescribed by A Moene for gabls1 (in dales)
        const double ce1 = 0.19;
        const double ce2 = 0.51;
        const double cm  = 0.12;
        const double ch1 = 1.0;
        const double ch2 = 2.0;
        const double cn  = 0.76;

        const double sgstkemin = 1e-9; // minimum value of SGS TKE to prevent model crash

        #ifdef USECUDA
        double* mlen_g;
        #endif
};
#endif
