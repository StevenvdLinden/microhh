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

#include <iostream>
#include "master.h"
#include "model.h"
#include "field3d.h"
#include "radiation_disabled.h"

Radiation_disabled::Radiation_disabled(Model* modelin, Input* inputin) : Radiation(modelin, inputin)
{
    swradiation = "0";
}

Radiation_disabled::~Radiation_disabled()
{
}

void Radiation_disabled::exec()
{
}

void Radiation_disabled::get_surface_radiation(Field3d* Qnet)
{
    master->print_error("Radiation_disabled can not provide surface radiation\n");
    throw 1;
}
