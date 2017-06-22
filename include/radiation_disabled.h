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

#ifndef RADIATON_DISABLED
#define RADIATON_DISABLED

#include "radiation.h"

class Model;
class Input;

/**
 * Derived class for a disabled radiation scheme
 */
class Radiation_disabled : public Radiation
{
    public:
        Radiation_disabled(Model*, Input*); ///< Constructor of the radiation class.
        ~Radiation_disabled();              ///< Destructor of the radiation class.

        void init();
        void create(Input*);
        void exec(); ///< Execute the radiation scheme.

        void get_surface_radiation(Field3d*);
};
#endif
