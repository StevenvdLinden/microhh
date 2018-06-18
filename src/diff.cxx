/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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
#include "grid.h"
#include "fields.h"
#include "master.h"
#include "defines.h"
#include "constants.h"
#include "model.h"
#include "stats.h"

// Moet hier ook al advection zichtbaar zijn? En automatisch ook derived class
//#include "advec.h"

// diffusion schemes
#include "diff.h"
#include "diff_disabled.h"
#include "diff_2.h"
#include "diff_4.h"
#include "diff_smag2.h"
#include "diff_scm.h"
#include "diff_sgs_tke.h"

Diff::Diff(Model* modelin, Input* inputin)
{
    model  = modelin;
    grid   = model->grid;
    fields = model->fields;
    master = model->master;

    swdiff = "0";

    int nerror = 0;
    nerror += inputin->get_item(&dnmax, "diff", "dnmax", "", 0.4);

    if (nerror)
        throw 1;
}

Diff::~Diff()
{
}

Diff* Diff::factory(Master* masterin, Input* inputin, Model* modelin, const std::string swspatialorder)
{
    std::string swdiff;
    std::string swboundary;
    std::string swadvec;

    int nerror = 0;
    nerror += inputin->get_item(&swdiff, "diff", "swdiff", "", swspatialorder);
    // load the boundary switch as well in order to be able to check whether the surface model is used
    nerror += inputin->get_item(&swboundary, "boundary", "swboundary", "", "default");
    // load the advection switch as well in order to be able to check whether it is disabled
    nerror += inputin->get_item(&swadvec, "advec", "swadvec", "", swspatialorder);

    if (nerror)
        return 0;

    if (swdiff == "0")
        return new Diff_disabled(modelin, inputin);
    else if (swdiff == "2")
        return new Diff_2(modelin, inputin);
    else if (swdiff == "4")
        return new Diff_4(modelin, inputin);
    else if (swdiff == "smag2")
    {
        // the subgrid model requires a surface model because of the MO matching at first level
        if ((swboundary == "surface") || (swboundary == "surface_bulk") || (swboundary == "surface_patch"))
            return new Diff_smag_2(modelin, inputin);
        else
        {
            masterin->print_error("swdiff=\"smag2\" requires a surface model (swboundary = \"surface\", \"surface_bulk\" or \"surface_patch\")\n");
            throw 1;
        }
    }
    else if (swdiff == "sgs_tke")
    {
        // the subgrid model requires a surface model because of the MO matching at first level
        if ((swboundary == "surface" && swspatialorder == "2")) //|| (swboundary == "surface_bulk") || (swboundary == "surface_patch"))
            return new Diff_sgs_tke(modelin, inputin);
        else
        {
            masterin->print_error("swdiff=\"sgs_tke\" currently only works with (swboundary = \"surface\") and/or with (swspatialorder = \"2\")\n");
            //masterin->print_error("swdiff=\"sgs_tke\" requires a surface model (swboundary = \"surface\", \"surface_bulk\" or \"surface_patch\")\n");
            throw 1;
        }
    }
    else if (swdiff == "scm")
    {

        masterin->print_error("swdiff=\"SCM\" is not implemented yet ;)\n");
        throw 1;

        // the single column model requires a surface model because of the MO matching at first level
        if ((swboundary == "surface") || (swboundary == "surface_bulk"))
            return new Diff_scm(modelin, inputin);
        else
        {
            masterin->print_error("swdiff=\"SCM\" requires a surface model (swboundary = \"surface\" or \"surface_bulk\")\n");
            throw 1;
        }

        if (!(swadvec == "0"))
        {
          masterin->print_error("swadvec=\"2\" or \"4\" is not allowed in single column mode; use swadvec=\"0\"\n");
          throw 1;
        }
    }
    else
    {
        masterin->print_error("\"%s\" is an illegal value for swdiff\n", swdiff.c_str());
        throw 1;
    }
}

std::string Diff::get_switch()
{
    return swdiff;
}
