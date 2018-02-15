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

#ifndef BOUNDARY_LSM
#define BOUNDARY_LSM

#include "boundary.h"
#include "stats.h"

class Model;
class Input;
class Stats;
struct Mask;

class Boundary_surface : public Boundary
{
    public:

        Boundary_surface(Model*, Input*);
        ~Boundary_surface();

        virtual void init(Input*);
        void create(Input*);
        virtual void set_values();

        void exec_stats(Mask*); ///< Execute statistics of surface
        void exec_cross();      ///< Execute cross sections of surface

        // Make these variables public for out-of-class usage.


#ifdef USECUDA
        // GPU functions and variables
        void prepare_device();
        void clear_device();
        void forward_device();  // TMP BVS
        void backward_device(); // TMP BVS
#endif

    protected:
        void process_input(Input *);   // Process and check the surface input
        void init_surface();           // Allocate and initialize the surface arrays
        void init_solver();            // Prepare the lookup table's for the surface layer solver
        void set_ustar();              // Set fixed ustar

    private:

        // surface scheme
        void update_bcs();///< Update the boundary values.
        void update_slave_bcs(); ///< Update the slave boundary values.

#ifdef USECUDA
        float* zL_sl_g;
        float* f_sl_g;
#endif
        int thermobc;

    protected:
        // cross sections
        std::vector<std::string> crosslist;        // List with all crosses from ini file
        std::vector<std::string> allowedcrossvars; // List with allowed cross variables

        Stats* stats;
        void update_slave_bcs();
};
#endif
