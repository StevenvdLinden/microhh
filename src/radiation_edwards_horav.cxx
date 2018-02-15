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

#include "radiation_edwards_horav.h"

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

Radiation_edwards_horav::Radiation_edwards_horav(Model* modelin, Input* inputin) : Radiation(modelin, inputin)
{
    swradiation = "edwards_horav";

    // copy the pointers
    grid    = model->grid;
    fields  = model->fields;
    stats   = model->stats;

    // set radiation pointers to zero
    lupf = 0;
    ldnf = 0;
    radt = 0;
    nofl = 0;

    // obtain specific input parameters
    int nerror = 0;
    nerror += inputin->get_item(&nlines, "radiation", "nlines", "");
    nerror += inputin->get_item(&raditer, "radiation", "raditer", "", 1000);
    nerror += inputin->get_item(&dtrad, "radiation", "dtrad", "", Constants::dbig);

    if (nerror)
        throw 1;

    // add field for radiation tendency
    //fields->init_diagnostic_field("radt", "Radiative tendency", "K s-1");

}

Radiation_edwards_horav::~Radiation_edwards_horav()
{
    delete[] lupf;
    delete[] ldnf;
    delete[] radt;
    delete[] Tabs;
    delete[] nofl;
}

void Radiation_edwards_horav::init(double ifactor)
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

    // container to (temporarily) store radiative fluxes
    lupf = new double[grid->kcells];
    ldnf = new double[grid->kcells];
    radt = new double[grid->kcells];
    Tabs = new double[grid->kcells];
    nofl = new double[grid->kcells];

    // specific integer sample time for radiation
    idtrad = (unsigned long)(ifactor * dtrad);
}

void Radiation_edwards_horav::create(Input* inputin)
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

unsigned long Radiation_edwards_horav::get_time_limit(unsigned long itime)
{
    if (swradiation == "0")
        return Constants::ulhuge;

    unsigned long idtlim = idtrad - itime % idtrad;

    return idtlim;
}

void Radiation_edwards_horav::exec()
{

    // check that if-statement below is correctly implemented, in combination with
    if (model->timeloop->get_itime() % idtrad == 0 )
    {

        // pointers to mean of temporary fields (at CPU)
        double* upflux   = fields->atmp["tmp1"]->datamean;
        double* dnflux   = fields->atmp["tmp1"]->datamean;
        double* bandflux = fields->atmp["tmp2"]->datamean;

        double* tempmean = fields->sp["th"]->datamean;
        grid->calc_mean(tempmean, fields->sp["th"]->data, grid->kcells);

        //model->thermo->get_thermo_field(fields->atmp["tmp1"], fields->atmp["tmp2"], "T_abs", true);
        //grid->calc_mean(Tabs, fields->atmp["tmp1"]->data, grid->kcells);

        for (int k=0; k<grid->kcells; ++k)
        {
            upflux[k] = 0;
        }

        // calculation of upward radiation fluxes and resulting tendencies
        if (grid->swspatialorder == "2")
        {
            calc_radiation_fluxes_up_2(upflux, bandflux, tempmean, grid->dzh);
        }
        else if (grid->swspatialorder == "4")
        {
            calc_radiation_fluxes_up_4(upflux, bandflux, tempmean, grid->dzh);
        }

        calc_radiation_tendency_up(radt, grid->dzh, upflux, lupf, arbc, crbc);

        for (int k=0; k<grid->kcells; ++k)
        {
            dnflux[k] = 0;
        }

        // calculation of downward radiation fluxes and resulting tendencies
        if (grid->swspatialorder == "2")
        {
            calc_radiation_fluxes_dn_2(dnflux, bandflux, tempmean, grid->dzh);
        }
        else if (grid->swspatialorder == "4")
        {
            calc_radiation_fluxes_dn_4(dnflux, bandflux, tempmean, grid->dzh);
        }

        calc_radiation_tendency_dn(radt, grid->dzh, dnflux, ldnf, arbc, crbc);

    }

    // apply the total radiative tendencies to the temperature tendencies.
    // 1) eventually radiation infrastructure to be used with all different thermo-modules:
    //    so a switch, call to proper temperature field needed here
    // 2) this function does not set the surface temperature tendency.
    //    Surface temperatures yet to be implemented via boundary-functions and Robin BC!

    apply_radiation_tendency(fields->st["th"]->data, radt);
}

void Radiation_edwards_horav::get_surface_radiation(Field3d* Qnet)
{
    master->print_error("Radiation_edwards can not provide surface radiation\n");
    throw 1;
}

void Radiation_edwards_horav::calc_radiation_fluxes_up_2(double* restrict flux_up, double* restrict fl_up_nr, const double* restrict temp, const double* restrict dzh) // totale upward and downward fluxen als argument meegeven, waar initialisatie?
{
    const int ii  = 1;
    const int jj  = grid->icells;
    const int kk  = grid->ijcells;
    const int ktot= grid->kcells;
    const int ntot= grid->ncells;

    const double tau_min = 1e-8; // constant minimal optical depth
    const double rho_air = 0.937; // air density for radiation
    const double diffus  = 1.66; // optical diffusivity

    double* Tbot = fields->sp["th"]->databot;

    // loop over radiation bands
    for (int kr=0; kr<nlines; ++kr)
    {
        // calculate outgoing radiative fluxe at surface
        int kb = grid->kstart;
        // DEZE IS WELLICHT NOG WEL FOUT!, werkt alleen omdat ik homogene bodemtemperatuur heb
        const int ij = grid->istart + grid->jstart*jj;

        // this has to be used for now, as Tabs at ghost cell would not be defined properly!
        fl_up_nr[kb] = pl0[kr] + pl1[kr] * Tbot[ij];
        // function above could be replaced by ??:
        // fl_up_nr[kb] = pl0[kr] + pl1[kr] * interp2(temp[kb],temp[kb-1]);

        // calculate upward fluxes
        for (int k=grid->kstart+1; k<(grid->kend+1); ++k)
        {
            double tau    = tau_min;    // (re)-initialize tau at tau_min
            double trans  = 0;          // (re)-initialize transmissivity
            double bbt    = 0, bbb = 0; // (re)-initialize planckian values

            // calculate optical depth and transmissivity
            tau  += ( c0[kr] + q_bg[k-1] * (c1[kr] + c2[kr] * q_bg[k-1]) ) * dzh[k-1] * rho_air;
            trans = std::exp(-diffus * tau);

            // planckian functions at top and bottom of air layers
            bbt = pl0[kr] + pl1[kr] * interp2(temp[k-1], temp[k]);
            bbb = pl0[kr] + pl1[kr] * interp2(temp[k-2], temp[k-1]);

            std::printf("%4d\t%7.4f\t%7.4f\n",k,interp2(temp[k-1], temp[k]),interp2(temp[k-2], temp[k-1]));

            // calculate upward flux for current level and radiation band
            fl_up_nr[k] = bbt + trans * (fl_up_nr[k-1] - bbb) - (bbt - bbb) * (1.0 - trans) / (diffus * tau);

        }

        // sum radiation bands : this instruction can be placed in the above loop. For clarity reasons, it is done here seperately
        for (int k=grid->kstart; k<(grid->kend+1); ++k)
        {
            flux_up[k] += fl_up_nr[k];
        }
    }
}

void Radiation_edwards_horav::calc_radiation_fluxes_dn_2(double* restrict flux_dn, double* restrict fl_dn_nr, const double* restrict temp, const double* restrict dzh) // totale upward and downward fluxen als argument meegeven, waar initialisatie?
{
    const int ii  = 1;
    const int jj  = grid->icells;
    const int kk  = grid->ijcells;
    const int ktot= grid->kcells;
    const int ntot= grid->ncells;

    const double tau_min = 1e-8; // constant minimal optical depth
    const double rho_air = 0.937; // air density for radiation
    const double diffus  = 1.66; // optical diffusivity

    // loop over radiation bands
    for (int kr=0; kr<nlines; ++kr)
    {
        // calculate incoming radiative fluxes at top
        int ke = grid->kend;
        fl_dn_nr[ke] = fdn_above_0[kr] + fdn_above_1[kr] * interp2(temp[ke],temp[ke-1]);

        // calculate downward fluxes
        for (int k=(grid->kend); k>(grid->kstart); --k) // Go explicitly to lowest level, situated at kstart
        {

            double tau    = tau_min;    // (re)-initialize tau at tau_min
            double trans  = 0;          // (re)-initialize transmissivity
            double bbt    = 0, bbb = 0; // (re)-initialize planckian values

            // calculate optical depth and transmissivity
            tau  += ( c0[kr] + q_bg[k-1] * (c1[kr] + c2[kr] * q_bg[k-1]) ) * dzh[k-1] * rho_air;
            trans = std::exp(-diffus * tau);

            // planckian functions at top and bottom of air layers
            bbt = pl0[kr] + pl1[kr] * interp2(temp[k-1], temp[k]);
            bbb = pl0[kr] + pl1[kr] * interp2(temp[k-2], temp[k-1]);

            // calculate downward flux for current level and radiation band
            fl_dn_nr[k-1] = bbb + trans * (fl_dn_nr[k] - bbt) - (bbb - bbt) * (1.0 - trans) / (diffus * tau) ; // -kk points to height lower?? CHECK grid counting !

        }

        // sum radiation bands
        for (int k=grid->kstart; k<(grid->kend+1); ++k)
        {
            flux_dn[k] += fl_dn_nr[k];
        }
    }
}

void Radiation_edwards_horav::calc_radiation_fluxes_up_4(double* restrict flux_up, double* restrict fl_up_nr, const double* restrict temp, const double* restrict dzh) // totale upward and downward fluxen als argument meegeven, waar initialisatie?
{
    const int ii  = 1;
    const int jj  = grid->icells;
    const int kk  = grid->ijcells;
    const int ntot= grid->ncells;

    const double tau_min = 1e-8;
    const double rho_air = 0.937;
    const double diffus  = 1.66;

    double * Tbot = fields->sp["th"]->databot;

    // loop over radiation bands
    for (int kr=0; kr<nlines; ++kr)
    {
        // calculate outgoing radiative fluxes at surface
        int kb = grid->kstart;
        // DEZE IS WELLICHT NOG WEL FOUT!, werkt alleen omdat ik homogene bodemtemperatuur heb
        const int ij = grid->istart + grid->jstart*jj;

        // this has to be used for now, as Tabs at ghost cell would not be defined properly!
        fl_up_nr[kb] = pl0[kr] + pl1[kr] * Tbot[ij];
        // function above could be replaced by ??:
        // fl_up_nr[kb] = pl0[kr] + pl1[kr] * interp4(temp[kb+1],temp[kb],temp[kb-1],temp[kb-2]);

        // calculate upward fluxes
        for (int k=grid->kstart+1; k<(grid->kend+1); ++k)
        {

            double tau   = tau_min;    // (re)-initialize tau at tau_min
            double trans = 0;          // (re)-initialize transmissivity
            double bbt   = 0, bbb = 0; // (re)-initialize planckian values

            // calculate new optical depth and transmissivity
            tau  += ( c0[kr] + q_bg[k-1] * (c1[kr] + c2[kr] * q_bg[k-1]) ) * dzh[k-1] * rho_air;
            trans = std::exp(-diffus * tau);

            // planckian functions at top and bottom of air layers
            bbt = pl0[kr] + pl1[kr] * interp4(temp[k-2], temp[k-1], temp[k], temp[k+1]);
            bbb = pl0[kr] + pl1[kr] * interp4(temp[k-3], temp[k-2], temp[k-1], temp[k]);

            // calculate flux for current level and radiation band
            fl_up_nr[k] = bbt + trans * (fl_up_nr[k-1] - bbb) - (bbt - bbb) * (1.0 - trans) / (diffus * tau);

        }

        // sum radiation bands
        for (int k=grid->kstart; k<(grid->kend+1); ++k)
        {
            flux_up[k] += fl_up_nr[k];
        }
    }
}

void Radiation_edwards_horav::calc_radiation_fluxes_dn_4(double* restrict flux_dn, double* restrict fl_dn_nr, const double* restrict temp, const double* restrict dzh) // totale upward and downward fluxen als argument meegeven, waar initialisatie?
{
    const int ii  = 1;
    const int jj  = grid->icells;
    const int kk  = grid->ijcells;
    const int ntot= grid->ncells;

    const double tau_min = 1e-8; // constant minimal optical depth
    const double rho_air = 0.937; // air density for radiation
    const double diffus  = 1.66; // optical diffusivity

    // loop over radiation bands
    for (int kr=0; kr<nlines; ++kr)
    {
        // calculate incoming radiative fluxes at top
        int ke = grid->kend;
        fl_dn_nr[ke] = fdn_above_0[kr] + fdn_above_1[kr] * interp4(temp[ke-2],temp[ke-1],temp[ke],temp[ke+1]);

        // calculate downward fluxes
        for (int k=(grid->kend); k>(grid->kstart); --k)
        {

            double tau   = tau_min;    // (re)-initialize tau at tau_min
            double trans = 0;          // (re)-initialize transmissivity
            double bbt   = 0, bbb = 0; // (re)-initialize planckian values

            // calculate new optical depth and transmissivity
            tau  += ( c0[kr] + q_bg[k-1] * (c1[kr] + c2[kr] * q_bg[k-1]) ) * dzh[k-1] * rho_air;
            trans = std::exp(-diffus * tau);

            // planckian functions at top and bottom of air layers
            bbt = pl0[kr] + pl1[kr] * interp4(temp[k-2], temp[k-1], temp[k], temp[k+1]);
            bbb = pl0[kr] + pl1[kr] * interp4(temp[k-3], temp[k-2], temp[k-1], temp[k]);

            // calculate flux for current level and radiation band
            fl_dn_nr[k-1] = bbb + trans * (fl_dn_nr[k] - bbt) - (bbb - bbt) * (1.0 - trans) / (diffus * tau) ; // -kk points to height lower?? CHECK grid counting !

        }

        // sum radiation bands
        for (int k=grid->kstart; k<(grid->kend+1); ++k)
        {
            flux_dn[k] += fl_dn_nr[k];
        }
    }
}

void Radiation_edwards_horav::calc_radiation_tendency_up(double* restrict radt, const double* restrict dzh, const double* restrict flux_up, double* restrict lupf, double* restrict QnT, double* restrict QnC)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const double rho_air = 0.937; // air density for radiation

    for (int k=grid->kstart; k<grid->kend; ++k)
    {

        // Quick fix to write longwave fluxes to statistics files
        lupf[k] = flux_up[k];
        radt[k] = - (flux_up[k+1] - flux_up[k]) / (Constants::cp * rho_air * dzh[k]);

    }

    // Quick fix to write longwave fluxes to statistics files
    lupf[grid->kend] = flux_up[grid->kend];

    // still do: these parameters have to be changed to set the Robin BC
    QnT = 0;
    QnC = 0;

}

void Radiation_edwards_horav::calc_radiation_tendency_dn(double* restrict radt, const double* restrict dzh, const double* restrict flux_dn, double* restrict ldnf, double* restrict QnT, double* restrict QnC)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const double rho_air = 0.937; // air density for radiation

    for (int k=grid->kstart; k<grid->kend; ++k)
    {

        // Quick fix to write longwave fluxes to statistics files
        ldnf[k] = flux_dn[k];

        // for efficiency: calculation of radiative tendencies is currently split for upward and downward fluxes.
        // this function, therefore, increments (!) the radiative tendency calculated in the function above.
        // also pay special attention to plus & minus signs, switched with respect to function above (source and sink) of energy for current layer
        radt[k] += ( flux_dn[k+1] - flux_dn[k]) / (Constants::cp * rho_air * dzh[k]);

    }

    // Quick fix to write longwave fluxes to statistics files
    ldnf[grid->kend] = flux_dn[grid->kend];

    // still do: these parameters have to be changed to set the Robin BC
    QnT += 0;
    QnC += 0;

}

void Radiation_edwards_horav::apply_radiation_tendency(double* restrict tempt, double *restrict radt)
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

                tempt[ijk] += radt[k];
            }
}

void Radiation_edwards_horav::init_stat()
{

    if (model->stats->get_switch() == "1")
    {
        model->stats->add_prof("radt", "Radiative tendency", "K s-1", "z"); ///< Radiative tendencies
        model->stats->add_prof("lupf", "Longwave upward flux", "W m-2", "zh"); ///< Radiation fluxes are located at vertical cell interfaces
        model->stats->add_prof("ldnf", "Longwave downward flux", "W m-2", "zh"); ///< Radiation fluxes are located at vertical cell interfaces
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


void Radiation_edwards_horav::exec_stats(Mask *m)
{
    const double NoOffset = 0.;

    // define the location
    const int sloc[] = {0,0,0};

    // calculate the mean
    //model->stats->calc_mean(m->profs["radt"].data, fields->sd["radt"]->data, NoOffset, sloc,
            //fields->atmp["tmp3"]->data, model->stats->nmask);

    model->stats->write_profile (radt, m->profs["radt"].data, model->stats->nmask);
    model->stats->write_profileh(lupf, m->profs["lupf"].data, model->stats->nmask);
    model->stats->write_profileh(ldnf, m->profs["ldnf"].data, model->stats->nmask);

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
