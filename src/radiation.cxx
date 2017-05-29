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

#include <cstdio>
#include <cmath>
#include <algorithm>
#include "master.h"
#include "model.h"

#include "radiation.h"
#include "radiation_disabled.h"

Radiation::Radiation(Model* modelin, Input* inputin)
{
    model  = modelin;
    master = model->master;

    swradiation = "0";
}

Radiation::~Radiation()
{
}

Radiation* Radiation::factory(Master* masterin, Input* inputin, Model* modelin)
{
    std::string swradiation;
    if (inputin->get_item(&swradiation, "radiation", "swradiation", "", "0"))
        throw 1;

    if (swradiation == "0")
        return new Radiation_disabled(modelin, inputin);
    else
    {
        masterin->print_error("\"%s\" is an illegal value for swradiation\n", swradiation.c_str());
        throw 1;
    }
}

std::string Radiation::get_switch()
{
    return swradiation;
}
