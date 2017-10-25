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
    double dtrad; // sample time for radiation : best placement here, or in header files?
    unsigned long idtrad; // integer sample time for radiation : best placement here, or in header files?

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
    nerror += inputin->get_item(&raditer, "radiation", "raditer", "", 1000);
    nerror += inputin->get_item(&dtrad, "radiation", "dtrad", "", Constants::dbig);

    if (nerror)
        throw 1;

    // add field for radiation tendency
    fields->init_diagnostic_field("radt", "Radiative tendency", "K s-1");

}

Radiation_edwards::~Radiation_edwards()
{
}

void Radiation_edwards::init(double ifactor)
{

    // container for background absolute humidity
    q_bg.resize(grid->kcells);

    // container for linearized radiation coefficients
    c0.resize(nlines);
    c1.resize(nlines);
    c2.resize(nlines);
    pl0.resize(nlines);
    pl1.resize(nlines);
    fdn_above_0.resize(nlines);
    fdn_above_1.resize(nlines);

    // specific integer sample time for radiation
    idtrad = (unsigned long)(ifactor * dtrad);
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

    // initialize 1D-profile statistics
    init_stat();

}

unsigned long Radiation_edwards::get_time_limit(unsigned long itime)
{
    if (swradiation == "0")
        return Constants::ulhuge;

    unsigned long idtlim = idtrad - itime % idtrad;

    return idtlim;
}

void Radiation_edwards::exec()
{

    // check that if-statement below is correctly implemented, in combination with
    if (model->timeloop->get_itime() % idtrad == 0 )
    {
        // pointers to temporary fields
        double* upflux   = fields->atmp["tmp1"]->data;
        double* dnflux   = fields->atmp["tmp1"]->data;
        double* bandflux = fields->atmp["tmp2"]->data;

        // explicitly initialize to zero!
        for (int n=0; n<grid->ncells; ++n)
        {
            upflux[n] = 0;
        }

        // calculation of upward radiation fluxes and resulting tendencies
        if (grid->swspatialorder == "2")
        {
            calc_radiation_fluxes_up_2(upflux, bandflux, fields->sp["th"]->data, grid->dzh);
        }
        else if (grid->swspatialorder == "4")
        {
            calc_radiation_fluxes_up_4(upflux, bandflux, fields->sp["th"]->data, grid->dzh);
        }

        calc_radiation_tendency_up(fields->sd["radt"]->data, grid->dzh, upflux, arbc, crbc);

        // explicitly initialize to zero!
        for (int n=0; n<grid->ncells; ++n)
        {
            dnflux[n] = 0;
        }

        // calculation of downward radiation fluxes and resulting tendencies
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

    // apply the total radiative tendencies to the temperature tendencies.
    // 1) Eventually radiation infrastructure to be used with all different thermo-modules:
    //    so a switch, call to proper temperature field needed here
    // 2) This function does not set the surface temperature tendency.
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

    const double tau_min = 1e-8; // constant minimal optical depth
    const double rho_air = 1.15; // air density for radiation
    const double diffus  = 1.66; // optical diffusivity

    double* Tbot = fields->sp["th"]->databot;

    // loop over radiation bands
    for (int kr=0; kr<nlines; ++kr)
    {
        // calculate outgoing radiative fluxe at surface
        int kb = grid->kstart;

        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk_bot = i + j*jj + kb*kk;
                const int ij      = i + j*jj;

                fl_up_nr[ijk_bot] = pl0[kr] + pl1[kr] * Tbot[ij];
                // function above could be replaced by :
                // fl_up_nr[ijk_bot] = pl0[kr] + pl1[kr] * interp2(temp[ijk_bot],temp[ijk_bot-kk]); ?
            }

        // calculate upward fluxes
        for (int k=grid->kstart+1; k<(grid->kend+1); ++k)
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ij  = i + j*jj;

                    double tau    = tau_min;    // (re)-initialize tau at tau_min
                    double trans  = 0;          // (re)-initialize transmissivity
                    double bbt    = 0, bbb = 0; // (re)-initialize planckian values

                    // calculate optical depth and transmissivity
                    tau  += ( c0[kr] + q_bg[k-1] * (c1[kr] + c2[kr] * q_bg[k-1]) ) * dzh[k-1] * rho_air;
                    trans = std::exp(-diffus * tau);

                    // planckian functions at top and bottom of air layers
                    bbt = pl0[kr] + pl1[kr] * interp2(temp[ijk-kk], temp[ijk]);
                    bbb = pl0[kr] + pl1[kr] * interp2(temp[ijk-2*kk], temp[ijk-kk]);

                    // calculate upward flux for current level and radiation band
                    fl_up_nr[ijk] = bbt + trans * (fl_up_nr[ijk-kk] - bbb) - (bbt - bbb) * (1.0 - trans) / (diffus * tau);

                }

        // sum radiation bands
        for (int k=grid->kstart; k<(grid->kend+1); ++k)
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    flux_up[ijk] += fl_up_nr[ijk];
                }
    }
}

void Radiation_edwards::calc_radiation_fluxes_dn_2(double* restrict flux_dn, double* restrict fl_dn_nr, const double* restrict temp, const double* restrict dzh) // totale upward and downward fluxen als argument meegeven, waar initialisatie?
{
    const int ii  = 1;
    const int jj  = grid->icells;
    const int kk  = grid->ijcells;
    const int ntot= grid->ncells;

    const double tau_min = 1e-8; // constant minimal optical depth
    const double rho_air = 1.15; // air density for radiation
    const double diffus  = 1.66; // optical diffusivity

    // loop over radiation bands
    for (int kr=0; kr<nlines; ++kr)
    {
        // calculate incoming radiative fluxes at top
        int ke = grid->kend;

        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk_top = i + j*jj + ke*kk;
                const int ij      = i + j*jj;

                fl_dn_nr[ijk_top] = fdn_above_0[kr] + fdn_above_1[kr] * interp2(temp[ijk_top],temp[ijk_top-kk]);
            }

        // calculate downward fluxes
        for (int k=(grid->kend); k>(grid->kstart); --k) // Go explicitly to lowest level, situated at kstart
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ij  = i + j*jj;

                    double tau    = tau_min;    // (re)-initialize tau at tau_min
                    double trans  = 0;          // (re)-initialize transmissivity
                    double bbt    = 0, bbb = 0; // (re)-initialize planckian values

                    // calculate optical depth and transmissivity
                    tau  += ( c0[kr] + q_bg[k-1] * (c1[kr] + c2[kr] * q_bg[k-1]) ) * dzh[k-1] * rho_air;
                    trans = std::exp(-diffus * tau);

                    // planckian functions at top and bottom of air layers
                    bbt = pl0[kr] + pl1[kr] * interp2(temp[ijk-kk], temp[ijk]);
                    bbb = pl0[kr] + pl1[kr] * interp2(temp[ijk-2*kk], temp[ijk-kk]);

                    // calculate downward flux for current level and radiation band
                    fl_dn_nr[ijk-kk] = bbb + trans * (fl_dn_nr[ijk] - bbt) - (bbb - bbt) * (1.0 - trans) / (diffus * tau) ; // -kk points to height lower?? CHECK grid counting !

                }

        // sum radiation bands
        for (int k=grid->kstart; k<(grid->kend+1); ++k)
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    flux_dn[ijk] += fl_dn_nr[ijk];
                }
    }
}

void Radiation_edwards::calc_radiation_fluxes_up_4(double* restrict flux_up, double* restrict fl_up_nr, const double* restrict temp, const double* restrict dzh) // totale upward and downward fluxen als argument meegeven, waar initialisatie?
{
    const int ii  = 1;
    const int jj  = grid->icells;
    const int kk  = grid->ijcells;
    const int ntot= grid->ncells;

    const double tau_min = 1e-8;
    const double rho_air = 1.15;
    const double diffus  = 1.66;

    double * Tbot = fields->sp["th"]->databot;

    // loop over radiation bands
    for (int kr=0; kr<nlines; ++kr)
    {
        // calculate outgoing radiative fluxes at surface
        int kb = grid->kstart;

        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk_bot = i + j*jj + kb*kk;
                const int ij      = i + j*jj;

                fl_up_nr[ijk_bot] = pl0[kr] + pl1[kr] * Tbot[ij];
                // function above could be replaced by :
                // fl_up_nr[ijk_bot] = pl0[kr] + pl1[kr] * interp4(temp[ijk_bot+kk],temp[ijk_bot],temp[ijk_bot-kk],temp[ijk_bot-2*kk]); ?
            }

        // calculate upward fluxes
        for (int k=grid->kstart+1; k<(grid->kend+1); ++k)
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ij  = i + j*jj;

                    double tau   = tau_min;    // (re)-initialize tau at tau_min
                    double trans = 0;          // (re)-initialize transmissivity
                    double bbt   = 0, bbb = 0; // (re)-initialize planckian values

                    // calculate new optical depth and transmissivity
                    tau  += ( c0[kr] + q_bg[k-1] * (c1[kr] + c2[kr] * q_bg[k-1]) ) * dzh[k-1] * rho_air;
                    trans = std::exp(-diffus * tau);

                    // planckian functions at top and bottom of air layers
                    bbt = pl0[kr] + pl1[kr] * interp4(temp[ijk-2*kk], temp[ijk-kk], temp[ijk], temp[ijk+kk]);
                    bbb = pl0[kr] + pl1[kr] * interp4(temp[ijk-3*kk], temp[ijk-2*kk], temp[ijk-kk], temp[ijk]);

                    // calculate flux for current level and radiation band
                    fl_up_nr[ijk] = bbt + trans * (fl_up_nr[ijk-kk] - bbb) - (bbt - bbb) * (1.0 - trans) / (diffus * tau);

                }

        // sum radiation bands
        for (int k=grid->kstart; k<(grid->kend+1); ++k)
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    flux_up[ijk] += fl_up_nr[ijk];
                }
    }
}

void Radiation_edwards::calc_radiation_fluxes_dn_4(double* restrict flux_dn, double* restrict fl_dn_nr, const double* restrict temp, const double* restrict dzh) // totale upward and downward fluxen als argument meegeven, waar initialisatie?
{
    const int ii  = 1;
    const int jj  = grid->icells;
    const int kk  = grid->ijcells;
    const int ntot= grid->ncells;

    const double tau_min = 1e-8; // constant minimal optical depth
    const double rho_air = 1.15; // air density for radiation
    const double diffus  = 1.66; // optical diffusivity

    // loop over radiation bands
    for (int kr=0; kr<nlines; ++kr)
    {
        // calculate incoming radiative fluxes at top
        int ke = grid->kend;

        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {

                const int ijk_top = i + j*jj + ke*kk;
                const int ij      = i + j*jj;

                fl_dn_nr[ijk_top] = fdn_above_0[kr] + fdn_above_1[kr] * interp4(temp[ijk_top-2*kk],temp[ijk_top-kk],temp[ijk_top],temp[ijk_top+kk]);
            }

        // calculate downward fluxes
        for (int k=(grid->kend); k>(grid->kstart); --k)
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ij  = i + j*jj;

                    double tau   = tau_min;    // (re)-initialize tau at tau_min
                    double trans = 0;          // (re)-initialize transmissivity
                    double bbt   = 0, bbb = 0; // (re)-initialize planckian values

                    // calculate new optical depth and transmissivity
                    tau  += ( c0[kr] + q_bg[k-1] * (c1[kr] + c2[kr] * q_bg[k-1]) ) * dzh[k-1] * rho_air;
                    trans = std::exp(-diffus * tau);

                    // planckian functions at top and bottom of air layers
                    bbt = pl0[kr] + pl1[kr] * interp4(temp[ijk-2*kk], temp[ijk-kk], temp[ijk], temp[ijk+kk]);
                    bbb = pl0[kr] + pl1[kr] * interp4(temp[ijk-3*kk], temp[ijk-2*kk], temp[ijk-kk], temp[ijk]);

                    // calculate flux for current level and radiation band
                    fl_dn_nr[ijk-kk] = bbb + trans * (fl_dn_nr[ijk] - bbt) - (bbb - bbt) * (1.0 - trans) / (diffus * tau) ; // -kk points to height lower?? CHECK grid counting !

                }

        // sum radiation bands
        for (int k=grid->kstart; k<(grid->kend+1); ++k)
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    flux_dn[ijk] += fl_dn_nr[ijk];
                }
    }
}

void Radiation_edwards::calc_radiation_tendency_up(double* restrict radtend, const double* restrict dzh, const double* restrict flux_up, double* restrict QnT, double* restrict QnC)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const double rho_air = 1.15; // air density for radiation

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // for efficiency: calculation of radiative tendencies is currently split for upward and downward fluxes.
                // this function overwrites (!) the radiative tendency from previous timesteps.
                radtend[ijk] = - (flux_up[ijk+kk] - flux_up[ijk]) / (Constants::cp * rho_air * dzh[k]);

            }

    // still do: these parameters have to be changed to set the Robin BC
    QnT = 0;
    QnC = 0;

}

void Radiation_edwards::calc_radiation_tendency_dn(double* restrict radtend, const double* restrict dzh, const double* restrict flux_dn, double* restrict QnT, double* restrict QnC)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const double rho_air = 1.15; // air density for radiation

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // for efficiency: calculation of radiative tendencies is currently split for upward and downward fluxes.
                // this function, therefore, increments (!) the radiative tendency calculated in the function above.
                // also pay special attention to plus & minus signs, switched with respect to function above (source and sink) of energy for current layer
                radtend[ijk] += ( flux_dn[ijk+kk] - flux_dn[ijk]) / (Constants::cp * rho_air * dzh[k]);

            }

    // still do: these parameters have to be changed to set the Robin BC
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
