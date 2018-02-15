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

#ifndef RADIATON_EDWARDS
#define RADIATON_EDWARDS

#include "radiation.h"

class Model;
class Input;
class Fields;
class Stats;
class Grid;

/**
 * Derived class for a disabled radiation scheme
 */
class Radiation_edwards: public Radiation
{
    public:
        Radiation_edwards(Model*, Input*); ///< Constructor of the radiation class.
        virtual ~Radiation_edwards();      ///< Destructor of the radiation class.

        void init(double);
        void create(Input*);
        void exec(); ///< Execute the radiation scheme.
        unsigned long get_time_limit(unsigned long); ///< Get the limit on the time step imposed by the radiation scheme.

        void get_surface_radiation(Field3d*);

        // Empty functions that are allowed to pass.
        void get_mask(Field3d*, Field3d*, Mask*) {}

        void exec_stats(Mask*);

    private:
        // Move to boundary_lsm.h later?
        double* arbc;
        double* crbc;

        void init_stat();

        void calc_radiation_fluxes_up_2(double *, double *, const double *, const double *);
        void calc_radiation_fluxes_dn_2(double *, double *, const double *, const double *);
        void calc_radiation_fluxes_up_4(double *, double *, const double *, const double *);
        void calc_radiation_fluxes_dn_4(double *, double *, const double *, const double *);
        void calc_intermediate_surf_rad(double *, const double *, double *, double *);
        void calc_radiation_tendency_up(double *, const double *, const double *, double *, double *);
        void calc_radiation_tendency_dn(double *, const double *, const double *, double *, double *);
        void apply_radiation_tendency(double *, double *);

        Stats* stats;

};
#endif
