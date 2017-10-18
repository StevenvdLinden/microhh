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

/*
 * ADDITIONAL NOTICE for radiation_edwards.cxx
 * This simplified radiation scheme is subject to the following copyright:
 *
 * British Crown Copyright (c) 2010, the Met Office
 * All rights reserved.
 * The use of this software is governed by the accompanying licence.
 * (nog ergens de licentie daadwerkelijk te plaatsen)
 */

#include <cstdio>
#include <algorithm>
#include <cmath>
#include <sstream>
#include "defines.h"
#include "master.h"
#include "model.h"
#include "field3d.h"
#include "constants.h"
#include "grid.h"
#include "fields.h"
#include "timeloop.h"
#include "stats.h"
#include "finite_difference.h"

#include "radiation_edwards.h"

namespace
{
    int nlines;
    int raditer;
    double dtrad;
    std::vector<double> q_bg;

    std::vector<double> c0;
    std::vector<double> c1;
    std::vector<double> c2;
    std::vector<double> pl0;
    std::vector<double> pl1;
    std::vector<double> fdn_above_0;
    std::vector<double> fdn_above_1;

}

using namespace Finite_difference::O2;
using namespace Finite_difference::O4;

Radiation_edwards::Radiation_edwards(Model* modelin, Input* inputin) : Radiation(modelin, inputin)
{
    swradiation = "edwards";

    // copy the pointers
    grid    = model->grid;
    fields  = model->fields;
    stats   = model->stats;

    // obtain specific input parameters
    int nerror = 0;
    nerror += inputin->get_item(&nlines, "radiation", "nlines", "");
    nerror += inputin->get_item(&raditer, "radiation", "raditer", "", 1000); // Default values high/low?
    nerror += inputin->get_item(&dtrad, "radiation", "dtrad", "", Constants::dbig); // Default values high/low?

    if (nerror)
        throw 1;

    // add field for radiation tendency
    fields->init_diagnostic_field("radt", "Radiative tendency", "K s-1"); // Check units
    // In deze setup is straling een diagnostisch veld, derhalve wordt deze dus niet opgeslagen voor restart times?
    // Geen probleem, maar het vereist dus altijd een herinitialisatie na een herstart!!

}

Radiation_edwards::~Radiation_edwards()
{
}

void Radiation_edwards::init()
{

    q_bg.resize(grid->kcells); // Container for background absolute humidity

    c0.resize(nlines);
    c1.resize(nlines);
    c2.resize(nlines);
    pl0.resize(nlines);
    pl1.resize(nlines);
    fdn_above_0.resize(nlines);
    fdn_above_1.resize(nlines);
}

void Radiation_edwards::create(Input* inputin)
{

    std::string data_file;

    int nerror = 0;
    const int kstart = grid->kstart;
    const int ktot = grid->ktot;
    const int ntot = grid->ncells;

    nerror += inputin->get_prof(&q_bg.data()[kstart], "q_bg", ktot);

    nerror += inputin->get_prof(c0.data(), "c0", nlines);
    nerror += inputin->get_prof(c1.data(), "c1", nlines);
    nerror += inputin->get_prof(c2.data(), "c2", nlines);
    nerror += inputin->get_prof(pl0.data(), "pl0", nlines);
    nerror += inputin->get_prof(pl1.data(), "pl1", nlines);
    nerror += inputin->get_prof(fdn_above_0.data(), "fdn_above_0", nlines);
    nerror += inputin->get_prof(fdn_above_1.data(), "fdn_above_1", nlines);

    if (nerror)
        throw 1;

    // initialiseer statistiek functions
    init_stat();

}

void Radiation_edwards::exec()
{
    int    iter;
    double time;

    // Get current iteration and time
    iter = model->timeloop->get_iteration();
    time = model->timeloop->get_time();

    // Let op: voor nu zijn arbc, en crbc gedeclareerd in radiation_edwards.h,
    // uiteindelijk moeten deze via boundary_lsm.h gedeclareerd worden en gecommuniceerd ?

    // Still add start of simulation, check for time
    // Dit moet slimmer voor het geval de timestep adaptief wordt ingesteld
    //if ((iter % raditer) == 0 )
    if (fmod(time, dtrad) == 0)
    {
        double* upflux   = fields->atmp["tmp1"]->data;
        double* dnflux   = fields->atmp["tmp1"]->data;
        double* bandflux = fields->atmp["tmp2"]->data;

        for (int n=0; n<grid->ncells; ++n)
        {
            upflux[n] = 0;
        }

        if (grid->swspatialorder == "2")
        {
            calc_radiation_fluxes_up_2(upflux, bandflux, fields->sp["th"]->data, grid->dzh);
        }
        else if (grid->swspatialorder == "4")
        {
            calc_radiation_fluxes_up_4(upflux, bandflux, fields->sp["th"]->data, grid->dzh);
        }

        calc_radiation_tendency_up(fields->sd["radt"]->data, grid->dzh, upflux, arbc, crbc);

        for (int n=0; n<grid->ncells; ++n)
        {
            dnflux[n] = 0;
        }

        if (grid->swspatialorder == "2")
        {
            calc_radiation_fluxes_dn_2(upflux, bandflux, fields->sp["th"]->data, grid->dzh);
        }
        else if (grid->swspatialorder == "4")
        {
            calc_radiation_fluxes_dn_4(upflux, bandflux, fields->sp["th"]->data, grid->dzh);
        }

        calc_radiation_tendency_dn(fields->sd["radt"]->data, grid->dzh, upflux, arbc, crbc);

    }

    // Apply the radiative tendencies to the (relevant) temperature tendencies.
    // 1) Eventually radiation infrastructure to be used with all different thermo-modules:
    //    so a switch, call to proper temperature field needed here
    // 2) This function does not set the surface temperature tendency: there is none.
    //    Surface temperatures yet to be implemented via boundary-functions and Robin BC!

    apply_radiation_tendency(fields->st["th"]->data,fields->sd["radt"]->data);
}

void Radiation_edwards::get_surface_radiation(Field3d* Qnet)
{
    master->print_error("Radiation_edwards can not provide surface radiation\n");
    throw 1;
}

void Radiation_edwards::calc_radiation_fluxes_up_2(double* restrict flux_up, double* restrict fl_up_nr, const double* restrict temp, const double* restrict dzh) // totale upward and downward fluxen als argument meegeven, waar initialisatie?
{
    const int ii  = 1;
    const int jj  = grid->icells;
    const int kk  = grid->ijcells;
    const int ntot= grid->ncells;

    const double tau_min = 1e-8; // Constant minimal optical depth
    const double rho_air = 1.15; // Air density for radiation --> specific per case, ask John? Or use diagnosed 'real' density?
    const double diffus = 1.66; // Optical diffusivity (what is the actual name?)

    double* Tbot = fields->sp["th"]->databot;

    for (int kr=0; kr<nlines; ++kr)
    {
        // Calculate radiative fluxes at boundaries
        int kb = grid->kstart;

        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk_bot = i + j*jj + kb*kk;
                const int ij      = i + j*jj;

                // 1] Nog de exner-functie doorgeven voor absolute temperatuur
                // 2] Mail v. John: nog niet geheel duidelijk voor mij hoe bijv. hogere emissiviteit oppervlak hierin is verwerkt

                // Outgoing radiation at the surface
                fl_up_nr[ijk_bot] = pl0[kr] + pl1[kr] * Tbot[ij];
                // Bovenstaande zou ook vervangen kunnen worden door? Klopt dit dan nog met updaten temperatuurrandvoorwaarde?:
                // fl_up_nr[ijk_bot] = pl0[kr] + pl1[kr] * interp2(temp[ijk_bot],temp[ijk_bot-kk]);
            }

        // 25/09/2017 : Er zijn flink wat indices veranderd om in lijn te komen met standaard. Ter verduidelijking:
        // upward flux tussen 1ste en 2de echte gridcell wordt opgeslagen op positie kstart+1, MAAR
        // vereist dus laagdikte en vochtgehalte van de laag eronder index = (k-1);
        // de temperatuur op te berekenen interface: interpolatie k en k-1;
        // en temperatuur op interface eronder: interpolatie k-1 en k-2 !

        // Calculate upward fluxes
        for (int k=grid->kstart+1; k<(grid->kend+1); ++k) // hier start bij kstart+1 : is dit correct?
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ij  = i + j*jj;

                    double tau = tau_min; // (Re)-Initialize tau, Overal moet toch tau_min bij opgeteld worden --> daarom in initialisatie
                    double trans = 0; // (Re)-Initialize transmissivity
                    double bbt = 0, bbb = 0; // (Re)-Initialize Planckian values

                    // Calculate new optical depth transmissivity
                    tau += ( c0[kr] + q_bg[k-1] * (c1[kr] + c2[kr] * q_bg[k-1]) ) * dzh[k-1] * rho_air;
                    trans = std::exp(-diffus * tau);

                    // Planckian functions at top and bottom of air layers
                    bbt = pl0[kr] + pl1[kr] * interp2(temp[ijk-kk], temp[ijk]);
                    bbb = pl0[kr] + pl1[kr] * interp2(temp[ijk-2*kk], temp[ijk-kk]);

                    // Increment total upward flux with upward flux for current radiation band
                    fl_up_nr[ijk] = bbt + trans * (fl_up_nr[ijk-kk] - bbb) - (bbt - bbb) * (1.0 - trans) / (diffus * tau);

                }

        // 25/09/2017 : Er zijn in feite 2 manieren om dit te doen, je laat de loop bij k=kend beginnen
        // en vervolgens bereken je de fluxen voor een niveau lager, of je begint bij kend-1, en berekent
        // de fluxen voor dat niveau gebruikmakend van de temperaturen op zowel niveau hoger als niveau lager.
        // Hier is voor de 1ste optie gekozen ook al loopt de loop dan wat anders,
        // omdat zo de structuur binnenin de loops gelijk zijn aan de berekening van de upward fluxen.

        // Sum upward fluxes
        for (int k=grid->kstart; k<(grid->kend+1); ++k) // Check, kstart=1, kend+1=34 --> 33 elements :P
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    flux_up[ijk] += fl_up_nr[ijk]; // Sum different radiation bands together
                }
    }

}

void Radiation_edwards::calc_radiation_fluxes_dn_2(double* restrict flux_dn, double* restrict fl_dn_nr, const double* restrict temp, const double* restrict dzh) // totale upward and downward fluxen als argument meegeven, waar initialisatie?
{
    const int ii  = 1;
    const int jj  = grid->icells;
    const int kk  = grid->ijcells;
    const int ntot= grid->ncells;

    const double tau_min = 1e-8; // Constant minimal optical depth
    const double rho_air = 1.15; // Air density for radiation --> specific per case, ask John? Or use diagnosed 'real' density?
    const double diffus = 1.66; // Optical diffusivity (what is the actual name?)

    for (int kr=0; kr<nlines; ++kr)
    {
        // Calculate radiative fluxes at boundaries
        int ke = grid->kend;

        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk_top = i + j*jj + ke*kk;
                const int ij      = i + j*jj;

                // 1] Nog de exner-functie doorgeven voor absolute temperatuur
                // 2] Mail v. John: nog niet geheel duidelijk voor mij hoe bijv. hogere emissiviteit oppervlak hierin is verwerkt

                // Incoming radiation at top op domain
                fl_dn_nr[ijk_top] = fdn_above_0[kr] + fdn_above_1[kr] * interp2(temp[ijk_top],temp[ijk_top-kk]);
            }

        // 25/09/2017 : Er zijn flink wat indices veranderd om in lijn te komen met standaard. Ter verduidelijking:
        // upward flux tussen 1ste en 2de echte gridcell wordt opgeslagen op positie kstart+1, MAAR
        // vereist dus laagdikte en vochtgehalte van de laag eronder index = (k-1);
        // de temperatuur op te berekenen interface: interpolatie k en k-1;
        // en temperatuur op interface eronder: interpolatie k-1 en k-2 !

        // 25/09/2017 : Er zijn in feite 2 manieren om dit te doen, je laat de loop bij k=kend beginnen
        // en vervolgens bereken je de fluxen voor een niveau lager, of je begint bij kend-1, en berekent
        // de fluxen voor dat niveau gebruikmakend van de temperaturen op zowel niveau hoger als niveau lager.
        // Hier is voor de 1ste optie gekozen ook al loopt de loop dan wat anders,
        // omdat zo de structuur binnenin de loops gelijk zijn aan de berekening van de upward fluxen.

        // Calculate downward fluxes
        for (int k=(grid->kend); k>(grid->kstart); --k) // Ga expliciet naar laagste level : deze zit dus op kstart
            for (int j=grid->jstart; j<grid->jend; ++j)
                //  #pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ij  = i + j*jj;

                    double tau = tau_min; // (Re)-Initialize tau, Overal moet toch tau_min bij opgeteld worden --> daarom in initialisatie
                    double trans = 0; // (Re)-Initialize transmissivity
                    double bbt = 0, bbb = 0; // (Re)-Initialize Planckian values

                    // Calculate new optical depth transmissivity
                    tau += ( c0[kr] + q_bg[k-1] * (c1[kr] + c2[kr] * q_bg[k-1]) ) * dzh[k-1] * rho_air;
                    trans = std::exp(-diffus * tau);

                    // Planckian functions at top and bottom of air layers
                    bbt = pl0[kr] + pl1[kr] * interp2(temp[ijk-kk], temp[ijk]);
                    bbb = pl0[kr] + pl1[kr] * interp2(temp[ijk-2*kk], temp[ijk-kk]);

                    // Increment total downward flux with downward flux for current radiation band
                    fl_dn_nr[ijk-kk] = bbb + trans * (fl_dn_nr[ijk] - bbt) - (bbb - bbt) * (1.0 - trans) / (diffus * tau) ; // -kk points to height lower?? CHECK grid counting !

                }

        // Sum downward fluxes
        for (int k=grid->kstart; k<(grid->kend+1); ++k) // Check, kstart=1, kend+1=34 --> 33 elementen :P
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    flux_dn[ijk] += fl_dn_nr[ijk]; // Sum different radiation bands together
                }
    }
}

// 25/09/2017 : Aan de 4de-orde implementatie wordt nog gewerkt

void Radiation_edwards::calc_radiation_fluxes_up_4(double* restrict flux_up, double* restrict fl_up_nr, const double* restrict temp, const double* restrict dzh) // totale upward and downward fluxen als argument meegeven, waar initialisatie?
{
    const int ii  = 1;
    const int jj  = grid->icells;
    const int kk  = grid->ijcells;
    const int ntot= grid->ncells;

    const double tau_min = 1e-8; // Constant minimal optical depth
    const double rho_air = 1.15; // Air density for radiation --> specific per case, ask John? Or use diagnosed 'real' density?
    const double diffus = 1.66; // Optical diffusivity (what is the actual name?)

    double * Tbot = fields->sp["th"]->databot;

    for (int kr=0; kr<nlines; ++kr)
    {
        // Calculate radiative fluxes at boundaries
        int kb = grid->kstart;
        int ke = grid->kend;

        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk_bot = i + j*jj + kb*kk;
                const int ij      = i + j*jj;

                // 1] Nog de exner-functie doorgeven voor absolute temperatuur
                // 2] Mail v. John: nog niet geheel duidelijk voor mij hoe bijv. hogere emissiviteit oppervlak hierin is verwerkt

                // Outgoing radiation at the surface
                fl_up_nr[ijk_bot] = pl0[kr] + pl1[kr] * Tbot[ij];
                // Bovenstaande zou ook vervangen kunnen worden door? Klopt dit dan nog met updaten temperatuurrandvoorwaarde?:
                // fl_up_nr[ijk_bot] = pl0[kr] + pl1[kr] * interp2(temp[ijk_bot],temp[ijk_bot-kk]);
            }

        // 25/09/2017 : Er zijn flink wat indices veranderd om in lijn te komen met standaard. Ter verduidelijking:
        // upward flux tussen 1ste en 2de echte gridcell wordt opgeslagen op positie kstart+1, MAAR
        // vereist dus laagdikte en vochtgehalte van de laag eronder index = (k-1);
        // de temperatuur op te berekenen interface: interpolatie k en k-1;
        // en temperatuur op interface eronder: interpolatie k-1 en k-2 !

        // Calculate upward fluxes
        for (int k=grid->kstart+1; k<(grid->kend+1); ++k) // hier start bij kstart+1 : is dit correct?
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ij  = i + j*jj;

                    double tau = tau_min; // (Re)-Initialize tau, Overal moet toch tau_min bij opgeteld worden --> daarom in initialisatie
                    double trans = 0; // (Re)-Initialize transmissivity
                    double bbt = 0, bbb = 0; // (Re)-Initialize Planckian values

                    // Calculate new optical depth transmissivity
                    tau += ( c0[kr] + q_bg[k-1] * (c1[kr] + c2[kr] * q_bg[k-1]) ) * dzh[k-1] * rho_air;
                    trans = std::exp(-diffus * tau);

                    // Planckian functions at top and bottom of air layers
                    bbt = pl0[kr] + pl1[kr] * interp4(temp[ijk-2*kk], temp[ijk-kk], temp[ijk], temp[ijk+kk]);
                    bbb = pl0[kr] + pl1[kr] * interp4(temp[ijk-3*kk], temp[ijk-2*kk], temp[ijk-kk], temp[ijk]);

                    // Increment total upward flux with upward flux for current radiation band
                    fl_up_nr[ijk] = bbt + trans * (fl_up_nr[ijk-kk] - bbb) - (bbt - bbb) * (1.0 - trans) / (diffus * tau);

                }

        // 25/09/2017 : Er zijn in feite 2 manieren om dit te doen, je laat de loop bij k=kend beginnen
        // en vervolgens bereken je de fluxen voor een niveau lager, of je begint bij kend-1, en berekent
        // de fluxen voor dat niveau gebruikmakend van de temperaturen op zowel niveau hoger als niveau lager.
        // Hier is voor de 1ste optie gekozen ook al loopt de loop dan wat anders,
        // omdat zo de structuur binnenin de loops gelijk zijn aan de berekening van de upward fluxen.

        // Sum upward fluxes
        for (int k=grid->kstart; k<(grid->kend+1); ++k) // Check, kstart=1, kend+1=34 --> 33 elementen :P
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    flux_up[ijk] += fl_up_nr[ijk]; // Sum different radiation bands together
                }
    }
}

void Radiation_edwards::calc_radiation_fluxes_dn_4(double* restrict flux_dn, double* restrict fl_dn_nr, const double* restrict temp, const double* restrict dzh) // totale upward and downward fluxen als argument meegeven, waar initialisatie?
{
    const int ii  = 1;
    const int jj  = grid->icells;
    const int kk  = grid->ijcells;
    const int ntot= grid->ncells;

    const double tau_min = 1e-8; // Constant minimal optical depth
    const double rho_air = 1.15; // Air density for radiation --> specific per case, ask John? Or use diagnosed 'real' density?
    const double diffus = 1.66; // Optical diffusivity (what is the actual name?)

    for (int kr=0; kr<nlines; ++kr)
    {
        // Calculate radiative fluxes at boundaries
        int ke = grid->kend;

        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {

                const int ijk_top = i + j*jj + ke*kk;
                const int ij      = i + j*jj;

                // 1] Nog de exner-functie doorgeven voor absolute temperatuur
                // 2] Mail v. John: nog niet geheel duidelijk voor mij hoe bijv. hogere emissiviteit oppervlak hierin is verwerkt

                // Incoming radiation at top op domain
                fl_dn_nr[ijk_top] = fdn_above_0[kr] + fdn_above_1[kr] * interp4(temp[ijk_top-2*kk],temp[ijk_top-kk],temp[ijk_top],temp[ijk_top+kk]);
            }

        // 25/09/2017 : Er zijn flink wat indices veranderd om in lijn te komen met standaard. Ter verduidelijking:
        // upward flux tussen 1ste en 2de echte gridcell wordt opgeslagen op positie kstart+1, MAAR
        // vereist dus laagdikte en vochtgehalte van de laag eronder index = (k-1);
        // de temperatuur op te berekenen interface: interpolatie k en k-1;
        // en temperatuur op interface eronder: interpolatie k-1 en k-2 !

        // 25/09/2017 : Er zijn in feite 2 manieren om dit te doen, je laat de loop bij k=kend beginnen
        // en vervolgens bereken je de fluxen voor een niveau lager, of je begint bij kend-1, en berekent
        // de fluxen voor dat niveau gebruikmakend van de temperaturen op zowel niveau hoger als niveau lager.
        // Hier is voor de 1ste optie gekozen ook al loopt de loop dan wat anders,
        // omdat zo de structuur binnenin de loops gelijk zijn aan de berekening van de upward fluxen.

        // Calculate downward fluxes
        for (int k=(grid->kend); k>(grid->kstart); --k) // Ga expliciet naar laagste level : deze zit dus op kstart
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ij  = i + j*jj;

                    double tau = tau_min; // (Re)-Initialize tau, Overal moet toch tau_min bij opgeteld worden --> daarom in initialisatie
                    double trans = 0; // (Re)-Initialize transmissivity
                    double bbt = 0, bbb = 0; // (Re)-Initialize Planckian values

                    // Calculate new optical depth transmissivity
                    tau += ( c0[kr] + q_bg[k-1] * (c1[kr] + c2[kr] * q_bg[k-1]) ) * dzh[k-1] * rho_air;
                    trans = std::exp(-diffus * tau);

                    // Planckian functions at top and bottom of air layers
                    bbt = pl0[kr] + pl1[kr] * interp4(temp[ijk-2*kk], temp[ijk-kk], temp[ijk], temp[ijk+kk]);
                    bbb = pl0[kr] + pl1[kr] * interp4(temp[ijk-3*kk], temp[ijk-2*kk], temp[ijk-kk], temp[ijk]);

                    // Increment total upward flux with upward flux for current radiation band
                    fl_dn_nr[ijk-kk] = bbb + trans * (fl_dn_nr[ijk] - bbt) - (bbb - bbt) * (1.0 - trans) / (diffus * tau) ; // -kk points to height lower?? CHECK grid counting !

                }

        // Sum downward fluxes
        for (int k=grid->kstart; k<(grid->kend+1); ++k) // Check, kstart=1, kend+1=34 --> 33 elementen :P
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    flux_dn[ijk] += fl_dn_nr[ijk]; // Sum different radiation bands together

                }
    }
}

void Radiation_edwards::calc_radiation_tendency_up(double* restrict radtend, const double* restrict dzh, const double* restrict flux_up, double* restrict QnT, double* restrict QnC)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const double rho_air = 1.15; // Air density for radiation --> specific per case, ask John? Or use diagnosed 'real' density?

    // Nog navragen/nadenken over hoe surface radiation mee te nemen. Deze moet meegegeven worden voor bepaling Robin-BC
    // Include surface radiation

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // Wat is eigenlijk de orde van deze berekening? Kan je over zoiets spreken bij straling
                // Ik neem aan dat straling alleen gebruikt kan worden met dry en/of moist thermodynamics,
                // als is mij nu niet geheel duidelijk wat dan anders is t.o.v. gebruik bouyancy parameter

                // For efficiency: calculation of radiative tendencies is currently split for upward and downward fluxes.
                // This function, therefore, overwrites (!) the radiative tendency from previous timesteps.

                radtend[ijk] = - (flux_up[ijk+kk] - flux_up[ijk]) / (Constants::cp * rho_air * dzh[k]); // Is dit juiste layer diepte? Is cp al gedefinieerd ergens?

            }

    QnT = 0;
    QnC = 0;

}

void Radiation_edwards::calc_radiation_tendency_dn(double* restrict radtend, const double* restrict dzh, const double* restrict flux_dn, double* restrict QnT, double* restrict QnC)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const double rho_air = 1.15; // Air density for radiation --> specific per case, ask John? Or use diagnosed 'real' density?

    // Nog navragen/nadenken over hoe surface radiation mee te nemen. Deze moet meegegeven worden voor bepaling Robin-BC
    // Include surface radiation

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // Wat is eigenlijk de orde van deze berekening? Kan je over zoiets spreken bij straling
                // Ik neem aan dat straling alleen gebruikt kan worden met dry en/of moist thermodynamics,
                // als is mij nu niet geheel duidelijk wat dan anders is t.o.v. gebruik bouyancy parameter

                // For efficiency: calculation of radiative tendencies is currently split for upward and downward fluxes.
                // This function, therefore, increments (!) the radiative tendency calculated in the function above.
                // Also pay special attention to plus & minus signs, switched with respect to function above (source and sink) of energy for current layer

                radtend[ijk] += ( flux_dn[ijk+kk] - flux_dn[ijk]) / (Constants::cp * rho_air * dzh[k]); // Is dit juiste layer diepte? Is cp al gedefinieerd ergens?

            }

    QnT += 0;
    QnC += 0;

}

void Radiation_edwards::apply_radiation_tendency(double* restrict tempt, double *restrict radt)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                tempt[ijk] += radt[ijk];
            }

}

void Radiation_edwards::init_stat()
{

    if (model->stats->get_switch() == "1")
    {
        model->stats->add_prof("radt", "Radiative tendency", "K s-1", "z");
        /* // Higher moments currently not needed?
           for (int n=2; n<5; ++n)
           {
           std::stringstream ss;
           ss << n;
           std::string sn = ss.str();
           stats->add_prof("radt"+sn, "Moment " +sn+" of the radiative tendence", "(K s-1)"+sn,"z");
           }
           */

    }
}


void Radiation_edwards::exec_stats(Mask *m)
{
    const double NoOffset = 0.;

    // define the location
    const int sloc[] = {0,0,0};

    // calculate the mean
    model->stats->calc_mean(m->profs["radt"].data, fields->sd["radt"]->data, NoOffset, sloc,
            fields->atmp["tmp3"]->data, model->stats->nmask);

    /*
    // calculate the moments
    for (int n=2; n<5; ++n)
    {
    std::stringstream ss;
    ss << n;
    std::string sn = ss.str();
    stats->calc_moment(fields->atmp["tmp1"]->data, m->profs["radt"].data, m->profs["radt"+sn].data, n, sloc,
    fields->atmp["tmp3"]->data, stats->nmask);
    }
    */

}

/*
   void Radiation_edwards::calc_radiation_fluxes_2(double* restrict temp, double* restrict dzh, double* restrict flux_up, double* restrict flux_dn) // totale upward and downward fluxen als argument meegeven, waar initialisatie?
   {
   const int ii  = 1;
   const int jj  = grid->icells;
   const int kk  = grid->ijcells;
   const int ntot= grid->ncells;

   const double tau_min = 1e-8; // Constant minimal optical depth
   const double rho_air = 1.15; // Air density for radiation --> specific per case, ask John? Or use diagnosed 'real' density?
   const double diffus = 1.66; // Optical diffusivity (what is the actual name?)

   double * Tbot = fields->sp["th"]->databot;

   for (int kr=0; kr<nlines; ++kr)
   {

   double fl_dn_nr[grid->ncells] ; // (Re)-Initialize flux fields for current radiation band
   double fl_up_nr[grid->ncells] ; // (Re)-Initialize flux fields for current radiation band

// Behandel surface en inkomend aan de top hier apart, kstart stelt surface voor, en kend absolute top??
int kb = grid->kstart-1;
int ke = grid->kend-1; // 17 Sept: laatste flux zou in

for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
for (int i=grid->istart; i<grid->iend; ++i)
{

// 29 augustus 2017: in deze indices zitten hoogstwaarschijnlijk de fouten, en worden de juist waardes dus niet aangeroepen!

const int ijk_bot = i + j*jj + kb*kk; // Is dit good coding practice ?
const int ijk_top = i + j*jj + ke*kk; // idem ?
const int ij      = i + j*jj;

// VRAGEN aan John: waarom is de straling aan het aardoppervlak ook afhankelijk van banden? En waarom geen standaard Boltzmann?
// Outgoing radiation at the surface
// HIER EVEN expliciet index ij geplaats!! 25 augustus 2017
// KLOPT Tbot[ij] hier? : is dit inderdaad de juiste temperatuur?
fl_up_nr[ijk_bot] = pl0[kr] + pl1[kr] * Tbot[ij]; //temp[ijk_bot]; // Insert absolute surface temperatuur hier

// Incoming radiation at top op domain
fl_dn_nr[ijk_top] = fdn_above_0[kr] + fdn_above_1[kr] * temp[ijk_top]; // Insert absolute temperature top of domain
}

// Calculate upward fluxes ; hier de resterende lagen
for (int k=grid->kstart; k<grid->kend; ++k) // REMOVED THE +1: kstart is first level ABOVE surface?? CHECK, 25 Aug 2017
for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
for (int i=grid->istart; i<grid->iend; ++i)
{
const int ijk = i + j*jj + k*kk;
const int ij  = i + j*jj;

double tau = tau_min; // (Re)-Initialize tau, Overal moet toch tau_min bij opgeteld worden --> daarom in initialisatie
double trans = 0; // (Re)-Initialize transmissivity
double bbt = 0, bbb = 0; // (Re)-Initialize Planckian values

// Calculate new optical depth transmissivity
tau += ( c0[kr] + q_bg[k] * (c1[kr] + c2[kr] * q_bg[k]) ) * dzh[k] * rho_air;
trans = std::exp(-diffus * tau);

if(k==grid->kstart)
{
// Planckian functions at top and bottom of air layers
bbt = pl0[kr] + pl1[kr] * temp[ijk]; //+ interp2(temp[ijk],temp[ijk+kk]) ofzoiets? // Insert absolute temperature here at top of cell
bbb = pl0[kr] + pl1[kr] * Tbot[ij]; // Wat is index voor 1 positie lager; kk toch? Insert absolute temperature here at bottom of cell
}
else
{
// Planckian functions at top and bottom of air layers
bbt = pl0[kr] + pl1[kr] * temp[ijk]; // Insert absolute temperature here at top of cell
bbb = pl0[kr] + pl1[kr] * temp[ijk-kk]; // Wat is index voor 1 positie lager; kk toch? Insert absolute temperature here at bottom of cell
}

// Increment total upward flux with upward flux for current radiation band
fl_up_nr[ijk] = bbt + trans * (fl_up_nr[ijk-kk] - bbb) - (bbt - bbb) * (1.0 - trans) / (diffus * tau) ; // -kk points to height lower?? CHECK grid counting !

}

// Calculate downward fluxes
for (int k=(grid->kend-1); k>(grid->kstart-1); --k) // GA EXPLICIET NAAR DIT LAAGSTE LEVEL // Move down through the layer, stop at last level above the surface == ksart+1??
for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
for (int i=grid->istart; i<grid->iend; ++i)
{

    const int ijk = i + j*jj + k*kk;
    const int ij  = i + j*jj;

    double tau = tau_min; // (Re)-Initialize tau, Overal moet toch tau_min bij opgeteld worden --> daarom in initialisatie
    double trans = 0; // (Re)-Initialize transmissivity
    double bbt = 0, bbb = 0; // (Re)-Initialize Planckian values

    // Calculate new optical depth transmissivity
    tau += ( c0[kr] + q_bg[k] * (c1[kr] + c2[kr] * q_bg[k]) ) * dzh[k] * rho_air;
    trans = std::exp(-diffus * tau);

    if(k==grid->kstart)
    {
        // Planckian functions at top and bottom of air layers
        bbt = pl0[kr] + pl1[kr] * temp[ijk]; // Insert absolute temperature here at top of cell
        bbb = pl0[kr] + pl1[kr] * Tbot[ij]; // Wat is index voor 1 positie lager; kk toch? Insert absolute temperature here at bottom of cell
    }
    else
    {
        // Planckian functions at top and bottom of air layers
        bbt = pl0[kr] + pl1[kr] * temp[ijk]; // Insert absolute temperature here at top of cell
        bbb = pl0[kr] + pl1[kr] * temp[ijk-kk]; // Wat is index voor 1 positie lager; kk toch? Insert absolute temperature here at bottom of cell
    }

    // Increment total upward flux with upward flux for current radiation band
    fl_dn_nr[ijk-kk] = bbb + trans * (fl_dn_nr[ijk] - bbt) - (bbb - bbt) * (1.0 - trans) / (diffus * tau) ; // -kk points to height lower?? CHECK grid counting !

}

// Sum upward and downward fluxes
for (int k=grid->kstart-1; k<grid->kend; ++k)
for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
for (int i=grid->istart; i<grid->iend; ++i)
{
    const int ijk = i + j*jj + k*kk;

    flux_up[ijk] += fl_up_nr[ijk]; // Sum different radiation bands together
    flux_dn[ijk] += fl_dn_nr[ijk]; // Sum different radiation bands together

}
}


// JUST FOR PRINTING FINAL ADDED RAD_FLUXES
for (int k=grid->kstart-1; k<grid->kend; ++k) // Note: counter starting at kstart (= suface?)
for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
for (int i=grid->istart; i<grid->iend; ++i)
{
    const int ijk = i + j*jj + k*kk;
    const int ijk_top = i + j*jj + (grid->kend)*kk; // idem ?
    std::printf("i=%d, j=%d, k=%d, ijk=%d ,rad_fluxes =\t%.14f\t%.14f\n",i,j,k,ijk,flux_up[ijk],flux_dn[ijk]);
}

}
*/
