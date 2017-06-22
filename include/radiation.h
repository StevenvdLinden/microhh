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

#ifndef RADIATON
#define RADIATON

#include <string>

class Master;
class Input;
class Model;
class Field3d;

/**
 * Base class for the radiation scheme. This class is abstract and only
 * derived classes can be instantiated. Derived classes are
 * implemented that handle different radiation schemes.
 */
class Radiation
{
    public:
        Radiation(Model*, Input*); ///< Constructor of the radiation class.
        virtual ~Radiation();      ///< Destructor of the radiation class.

        static Radiation* factory(Master*, Input*, Model*); ///< Factory function for radiation class generation.

        std::string get_switch();

        // Pure virtual functions that have to be implemented in derived class.
        virtual void init() = 0;   ///< Initialize the radiation scheme.
        virtual void create(Input*) = 0; ///< Create the data fields.
        virtual void exec() = 0;   ///< Execute the radiation scheme

        // Get en set interfaces
        virtual void get_surface_radiation(Field3d*) = 0;

    protected:
        Master* master; ///< Pointer to master class.
        Model*  model;  ///< Pointer to model class.

        std::string swradiation;
};
#endif
