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

#ifndef DIFF_SCM
#define DIFF_SCM

#include "diff.h"

class Diff_scm : public Diff
{
    public:
        Diff_scm(Model*, Input*);
        ~Diff_scm();

        void exec();
        void exec_viscosity();

        unsigned long get_time_limit(unsigned long, double);
        double get_dn(double);

        double tPr;

        #ifdef USECUDA
        // GPU functions and variables
        void prepare_device();
        void clear_device();
        #endif

        // Empty functions, there are allowed to pass.
        void set_values() {}

    private:
        template<bool>
        void calc_strain2(double*,
                          double*, double*, double*,
                          double*, double*,
                          double*, double*,
                          double*, double*, double*);

        void calc_evisc(double*,
                        double*, double*, double*, double*,
                        double*, double*, double*,
                        double*, double*,
                        double*, double*, double*,
                        double);

        template<bool>
        void calc_evisc_neutral(double*,
                                double*, double*, double*,
                                double*, double*,
                                double*, double*,
                                double, double);

        template<bool>
        void diff_u(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
        template<bool>
        void diff_v(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

        void diff_w(double*, double*, double*, double*, double*, double*, double*, double*, double*);
        void diff_c(double*, double*, double*, double*, double*, double*, double*, double*, double*, double);

        double calc_dnmul(double*, double*, double);

        double cs;

        #ifdef USECUDA
        double* mlen_g;
        #endif
};
#endif
