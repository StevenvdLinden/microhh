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

// NOOT: algemene vraag, waarom is op een aantal punten een functie voor resolved boundaries ingevoerd.
// in diff.cxx wordt een check gedaan op het gebruik van swboundary = surface

// NOOT : er zijn geen diffusie init/create-functies, dus nu is init_stat geplaatst in de set_values-functie (niet de netste oplossing)
// in alle andere diff-instances worden exec_stat() {} als empty function meegegeven

#include <cstdio>
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "master.h"
#include "diff_sgs_tke.h"
#include "boundary_surface.h"
#include "defines.h"
#include "constants.h"
#include "thermo.h"
#include "model.h"
#include "monin_obukhov.h"

#include "advec.h"

namespace
{
    namespace most = Monin_obukhov;
}

Diff_sgs_tke::Diff_sgs_tke(Model* modelin, Input* inputin) : Diff(modelin, inputin)
{
    swdiff = "sgs_tke";

    #ifdef USECUDA
    mlen_g = 0;
    #endif

    fields->init_prognostic_field("sgs_tke", "Subgrid scale TKE", "m2 s-2");
    fields->init_diagnostic_field("evisc", "Eddy viscosity", "m2 s-1");

    int nerror = 0;
    nerror += inputin->get_item(&dnmax, "diff", "dnmax", "", 0.5  );
    nerror += inputin->get_item(&cs   , "diff", "cs"   , "", 0.23 );

    if (nerror)
        throw 1;
}

Diff_sgs_tke::~Diff_sgs_tke()
{
#ifdef USECUDA
    clear_device();
#endif
}

void Diff_sgs_tke::set_values() // this is probably the ugliest solution, there is no create function for the diffusion class, so just initialise stats here ...
{
    if (model->thermo->get_switch() != "0")
        init_stats();
}

#ifndef USECUDA
unsigned long Diff_sgs_tke::get_time_limit(const unsigned long idt, const double dt)
{
    // Ugly solution for now, to avoid passing entire tPr-field, or recalculate it again; SvdLinden, May 2018
    double Pr_max = 2.0;

    double dnmul = calc_dnmul(fields->sd["evisc"]->data, grid->dzi, Pr_max);
    // Avoid zero division.
    dnmul = std::max(Constants::dsmall, dnmul);

    return idt * dnmax/(dnmul*dt);
}
#endif

#ifndef USECUDA
double Diff_sgs_tke::get_dn(const double dt)
{
    // Ugly solution for now, to avoid passing entire tPr-field, or recalculate it again; SvdLinden, May 2018
    double Pr_max = 2.0;

    // calculate eddy viscosity
    const double dnmul = calc_dnmul(fields->sd["evisc"]->data, grid->dzi, Pr_max);

    return dnmul*dt;
}
#endif

#ifndef USECUDA

// NOOT: voor start tijdsintegratie en na elke integratiestap wordt via boundary->exec() alle BCs bijgewerkt, i.e., worden de ghostcells bijgewerkt
// Zoek nog uit: ergens worden hier steeds de boundary_surface functies gesubstituted wanneer surface-model wordt gebruikt (via de factory)
// maar bijv. boundary_surface.cxx heeft geen exec()-functie: hier wordt dus de exec uit standaard boundary.cxx gebruikt??

// NOOT: in het algemeen worden er in een aantal functies op dit moment nutteloze parameters meegegeven, zoals boundaryptr, of datafluxbots wanneer wall resolved is?
// waarom niet overal nullpointers meegeven dan?

void Diff_sgs_tke::exec_viscosity()
{
    // Do a cast because the base boundary class does not have the MOST related variables.
    Boundary_surface* boundaryptr = static_cast<Boundary_surface*>(model->boundary);

    // NOOT : zover ik nu zie, is dit in zijn geheel niet nodig voor berekening van de viscositeiten

    // // Calculate strain rate using MO for velocity gradients lowest level
    // if (model->boundary->get_switch() == "surface")
    //     calc_strain2<false>(fields->sd["evisc"]->data,
    //                         fields->u->data, fields->v->data, fields->w->data,
    //                         fields->u->datafluxbot, fields->v->datafluxbot,
    //                         boundaryptr->ustar, boundaryptr->obuk,
    //                         grid->z, grid->dzi, grid->dzhi);
    // // Calculate strain rate using resolved boundaries
    // else
    //     calc_strain2<true>(fields->sd["evisc"]->data,
    //                        fields->u->data, fields->v->data, fields->w->data,
    //                        fields->u->datafluxbot, fields->v->datafluxbot,
    //                        NULL, NULL, // BvS, for now....
    //                        grid->z, grid->dzi, grid->dzhi);

    // start with retrieving the stability information
    if (model->thermo->get_switch() == "0")
    {
        // Calculate eddy viscosity using MO at lowest model level
        if (model->boundary->get_switch() == "surface")
            calc_evisc_neutral<false>(fields->sd["evisc"]->data,
                                      fields->sp["sge_tke"]->data, fields->u->data, fields->v->data, fields->w->data,
                                      fields->u->datafluxbot, fields->v->datafluxbot,
                                      grid->z, grid->dz, boundaryptr->z0m, fields->visc);
        // Calculate eddy viscosity assuming resolved walls
        else
            calc_evisc_neutral<true>(fields->sd["evisc"]->data,
                                     fields->sp["sge_tke"]->data, fields->u->data, fields->v->data, fields->w->data,
                                     fields->u->datafluxbot, fields->v->datafluxbot,
                                     grid->z, grid->dz, 0, fields->visc); // BvS, for now....
    }
    // assume buoyancy calculation is needed
    else
    {
        // NOOT: waarom wordt hier tmp2 voor tijdelijke berekeningen meegegeven? In de thermo-functie zelf wordt gewoon het huidige temperatuurveld meegegeven?
        // store the buoyancyflux in tmp1
        model->thermo->get_buoyancy_fluxbot(fields->atmp["tmp1"]);
        // retrieve the full field in tmp1 and use tmp2 for temporary calculations
        model->thermo->get_thermo_field(fields->atmp["tmp1"], fields->atmp["tmp2"], "N2", false);
        // model->thermo->getThermoField(fields->sd["tmp1"], fields->sd["tmp2"], "b");

        if (model->boundary->get_switch() == "surface")
            calc_evisc<false>(fields->sd["evisc"]->data,
                             fields->sp["sge_tke"]->data, fields->u->data, fields->v->data, fields->w->data, fields->atmp["tmp1"]->data,
                             fields->u->datafluxbot, fields->v->datafluxbot, fields->atmp["tmp1"]->datafluxbot,
                             boundaryptr->ustar, boundaryptr->obuk,
                             grid->z, grid->dz, grid->dzi,
                             boundaryptr->z0m, fields->visc);

        else
            calc_evisc<true>(fields->sd["evisc"]->data,
                            fields->sp["sge_tke"]->data, fields->u->data, fields->v->data, fields->w->data, fields->atmp["tmp1"]->data,
                            fields->u->datafluxbot, fields->v->datafluxbot, fields->atmp["tmp1"]->datafluxbot,
                            0, 0,
                            grid->z, grid->dz, grid->dzi,
                            0, fields->visc); // SJAvdLinden, for now just pass null pointers
    }
}
#endif

#ifndef USECUDA
void Diff_sgs_tke::exec()
{
    // Do a cast because the base boundary class does not have the MOST related variables.
    Boundary_surface* boundaryptr = static_cast<Boundary_surface*>(model->boundary);

    // pointers to temporary fields; SvdLinden, May 2018
    double* S2 = fields->atmp["tmp1"]->data;
    double* N2 = fields->atmp["tmp1"]->data;

    // Do calculation in this order so temp1 'can be recycled'; SvdLinden, May 2018

    // Calculate strain rate squared for use in SGS TKE production and store in tmp1
    if (model->boundary->get_switch() == "surface")
        calc_strain2<false>(S2,
                            fields->u->data, fields->v->data, fields->w->data,
                            fields->u->datafluxbot, fields->v->datafluxbot,
                            boundaryptr->ustar, boundaryptr->obuk,
                            grid->z, grid->dzi, grid->dzhi);
    // Calculate strain rate using resolved boundaries
    else
        calc_strain2<true>(S2,
                           fields->u->data, fields->v->data, fields->w->data,
                           fields->u->datafluxbot, fields->v->datafluxbot,
                           NULL, NULL, // BvS, for now....
                           grid->z, grid->dzi, grid->dzhi);

    calc_sgs_tke_shear_tend_2nd(fields->st["sgs_tke"]->data, fields->sp["evisc"]->data, S2);

    // Retrieve buoyancy field for use in SGS TKE production and store in tmp1
    // these functions have been exactly copied from the exec_viscosity()-function above
    // store the buoyancyflux in tmp1 (NOOT: WAAROM DEZE STAP?)
    model->thermo->get_buoyancy_fluxbot(fields->atmp["tmp1"]);
    // retrieve the full field in tmp1 and use tmp2 for temporary calculations
    model->thermo->get_thermo_field(fields->atmp["tmp1"], fields->atmp["tmp2"], "N2", false);

    calc_sgs_tke_buoyancy_tend_2nd(fields->st["sgs_tke"]->data, fields->sp["sgs_tke"]->data, fields->sp["evisc"]->data, N2,
                                  grid->dz);
    calc_sgs_tke_dissipation_2nd(fields->st["sgs_tke"]->data, fields->sp["sgs_tke"]->data, N2,
                                grid->dz);

    if(model->boundary->get_switch() == "surface")
    {
        diff_u<false>(fields->ut->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi, grid->dzhi, fields->sd["evisc"]->data,
               fields->u->datafluxbot, fields->u->datafluxtop, fields->rhoref, fields->rhorefh);
        diff_v<false>(fields->vt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi, grid->dzhi, fields->sd["evisc"]->data,
               fields->v->datafluxbot, fields->v->datafluxtop, fields->rhoref, fields->rhorefh);
        diff_w(fields->wt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi, grid->dzhi, fields->sd["evisc"]->data,
               fields->rhoref, fields->rhorefh);

        // SGS TKE field should diffuse with eddy viscosity for momentum, not for heat/scalars; SvdLinden, May 2018
        // NOOT: currently, this would work because tmp2 which holds N2 should NOT be overwritten already at this stage.
        for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); ++it)
        {
            if(it->second->name.c_str()== "sgs_tke") // NOOT : is dit de correcte verwijzing naar sgs_tke-veld, hoe werkt de Map ? Waarom ->second i.p.v. ->first ?
            {
                diff_sgs_tke<false>(it->second->data, fields->sp[it->first]->data, grid->dzi, grid->dzhi, fields->sd["evisc"]->data,
                             fields->sp[it->first]->datafluxbot, fields->sp[it->first]->datafluxtop, fields->rhoref, fields->rhorefh);
            }
            else
            {
                diff_c<false>(it->second->data, fields->sp[it->first]->data, grid->dz, grid->dzi, grid->dzhi, fields->sd["evisc"]->data, fields->sp["sgs_tke"]->data, N2,
                      fields->sp[it->first]->datafluxbot, fields->sp[it->first]->datafluxtop, fields->rhoref, fields->rhorefh);
            }
        }
    }
    else
    {
        diff_u<true>(fields->ut->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi, grid->dzhi, fields->sd["evisc"]->data,
               fields->u->datafluxbot, fields->u->datafluxtop, fields->rhoref, fields->rhorefh);
        diff_v<true>(fields->vt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi, grid->dzhi, fields->sd["evisc"]->data,
               fields->v->datafluxbot, fields->v->datafluxtop, fields->rhoref, fields->rhorefh);
        diff_w(fields->wt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi, grid->dzhi, fields->sd["evisc"]->data,
               fields->rhoref, fields->rhorefh);

        // SGS TKE field should diffuse with eddy viscosity for momentum, not for heat/scalars; SvdLinden, May 2018
        // NOOT: currently, this would work because tmp2 which holds N2 should NOT be overwritten already at this stage.
        for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); ++it)
        {
            if(it->second->name.c_str() == "sgs_tke")
            {
                diff_sgs_tke<true>(it->second->data, fields->sp[it->first]->data, grid->dzi, grid->dzhi, fields->sd["evisc"]->data,
                            fields->sp[it->first]->datafluxbot, fields->sp[it->first]->datafluxtop, fields->rhoref, fields->rhorefh);
            }
            else
            {
                diff_c<true>(it->second->data, fields->sp[it->first]->data, grid->dz, grid->dzi, grid->dzhi, fields->sd["evisc"]->data, fields->sp["sgs_tke"]->data, N2,
                      fields->sp[it->first]->datafluxbot, fields->sp[it->first]->datafluxtop, fields->rhoref, fields->rhorefh);
            }
        }
    }

}
#endif

template <bool resolved_wall>
void Diff_sgs_tke::calc_strain2(double* restrict strain2,
                               double* restrict u, double* restrict v, double* restrict w,
                               double* restrict ufluxbot, double* restrict vfluxbot,
                               double* restrict ustar, double* restrict obuk,
                               double* restrict z, double* restrict dzi, double* restrict dzhi)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    const int k_offset = resolved_wall ? 0 : 1;

    // If the wall isn't resolved, calculate du/dz and dv/dz at lowest grid height using MO
    if (!resolved_wall)
    {
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;

                strain2[ijk] = 2.*(
                               // du/dx + du/dx
                               + std::pow((u[ijk+ii]-u[ijk])*dxi, 2)

                               // dv/dy + dv/dy
                               + std::pow((v[ijk+jj]-v[ijk])*dyi, 2)

                               // dw/dz + dw/dz
                               + std::pow((w[ijk+kk]-w[ijk])*dzi[kstart], 2)

                               // du/dy + dv/dx
                               + 0.125*std::pow((u[ijk      ]-u[ijk   -jj])*dyi  + (v[ijk      ]-v[ijk-ii   ])*dxi, 2)
                               + 0.125*std::pow((u[ijk+ii   ]-u[ijk+ii-jj])*dyi  + (v[ijk+ii   ]-v[ijk      ])*dxi, 2)
                               + 0.125*std::pow((u[ijk   +jj]-u[ijk      ])*dyi  + (v[ijk   +jj]-v[ijk-ii+jj])*dxi, 2)
                               + 0.125*std::pow((u[ijk+ii+jj]-u[ijk+ii   ])*dyi  + (v[ijk+ii+jj]-v[ijk   +jj])*dxi, 2)

                               // du/dz
                               + 0.5*std::pow(-0.5*(ufluxbot[ij]+ufluxbot[ij+ii])/(Constants::kappa*z[kstart]*ustar[ij])*most::phim(z[kstart]/obuk[ij]), 2)

                               // dw/dx
                               + 0.125*std::pow((w[ijk      ]-w[ijk-ii   ])*dxi, 2)
                               + 0.125*std::pow((w[ijk+ii   ]-w[ijk      ])*dxi, 2)
                               + 0.125*std::pow((w[ijk   +kk]-w[ijk-ii+kk])*dxi, 2)
                               + 0.125*std::pow((w[ijk+ii+kk]-w[ijk   +kk])*dxi, 2)

                               // dv/dz
                               + 0.5*std::pow(-0.5*(vfluxbot[ij]+vfluxbot[ij+jj])/(Constants::kappa*z[kstart]*ustar[ij])*most::phim(z[kstart]/obuk[ij]), 2)

                               // dw/dy
                               + 0.125*std::pow((w[ijk      ]-w[ijk-jj   ])*dyi, 2)
                               + 0.125*std::pow((w[ijk+jj   ]-w[ijk      ])*dyi, 2)
                               + 0.125*std::pow((w[ijk   +kk]-w[ijk-jj+kk])*dyi, 2)
                               + 0.125*std::pow((w[ijk+jj+kk]-w[ijk   +kk])*dyi, 2) );

                // add a small number to avoid zero divisions
                strain2[ijk] += Constants::dsmall;
            }
    }

    for (int k=grid->kstart+k_offset; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                strain2[ijk] = 2.*(
                               // du/dx + du/dx
                               + std::pow((u[ijk+ii]-u[ijk])*dxi, 2)

                               // dv/dy + dv/dy
                               + std::pow((v[ijk+jj]-v[ijk])*dyi, 2)

                               // dw/dz + dw/dz
                               + std::pow((w[ijk+kk]-w[ijk])*dzi[k], 2)

                               // du/dy + dv/dx
                               + 0.125*std::pow((u[ijk      ]-u[ijk   -jj])*dyi  + (v[ijk      ]-v[ijk-ii   ])*dxi, 2)
                               + 0.125*std::pow((u[ijk+ii   ]-u[ijk+ii-jj])*dyi  + (v[ijk+ii   ]-v[ijk      ])*dxi, 2)
                               + 0.125*std::pow((u[ijk   +jj]-u[ijk      ])*dyi  + (v[ijk   +jj]-v[ijk-ii+jj])*dxi, 2)
                               + 0.125*std::pow((u[ijk+ii+jj]-u[ijk+ii   ])*dyi  + (v[ijk+ii+jj]-v[ijk   +jj])*dxi, 2)

                               // du/dz + dw/dx
                               + 0.125*std::pow((u[ijk      ]-u[ijk   -kk])*dzhi[k  ] + (w[ijk      ]-w[ijk-ii   ])*dxi, 2)
                               + 0.125*std::pow((u[ijk+ii   ]-u[ijk+ii-kk])*dzhi[k  ] + (w[ijk+ii   ]-w[ijk      ])*dxi, 2)
                               + 0.125*std::pow((u[ijk   +kk]-u[ijk      ])*dzhi[k+1] + (w[ijk   +kk]-w[ijk-ii+kk])*dxi, 2)
                               + 0.125*std::pow((u[ijk+ii+kk]-u[ijk+ii   ])*dzhi[k+1] + (w[ijk+ii+kk]-w[ijk   +kk])*dxi, 2)

                               // dv/dz + dw/dy
                               + 0.125*std::pow((v[ijk      ]-v[ijk   -kk])*dzhi[k  ] + (w[ijk      ]-w[ijk-jj   ])*dyi, 2)
                               + 0.125*std::pow((v[ijk+jj   ]-v[ijk+jj-kk])*dzhi[k  ] + (w[ijk+jj   ]-w[ijk      ])*dyi, 2)
                               + 0.125*std::pow((v[ijk   +kk]-v[ijk      ])*dzhi[k+1] + (w[ijk   +kk]-w[ijk-jj+kk])*dyi, 2)
                               + 0.125*std::pow((v[ijk+jj+kk]-v[ijk+jj   ])*dzhi[k+1] + (w[ijk+jj+kk]-w[ijk   +kk])*dyi, 2) );

                       // Add a small number to avoid zero divisions.
                       strain2[ijk] += Constants::dsmall;
            }

}

template <bool resolved_wall>
void Diff_sgs_tke::calc_evisc(double* restrict evisc,
                             double* restrict sgstke, double* restrict u, double* restrict v, double* restrict w,  double* restrict N2,
                             double* restrict ufluxbot, double* restrict vfluxbot, double* restrict bfluxbot,
                             double* restrict ustar, double* restrict obuk,
                             double* restrict z, double* restrict dz, double* restrict dzi,
                             const double z0m, const double mvisc)
{
    // Variables for the length scales.
    double mlen,mlen0,fac;

    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    const double dx = grid->dx;
    const double dy = grid->dy;

    // local copies to aid vectorization
    // NOOT: volgens mij is het maken van local copies onnodig voor consts ?
    double cm  = this->cm;
    double cn  = this->cn;

    // NOOT , anyway: waarom wordt hier lokaal MO gebruikt voor de RitPrratio, en 'mist' bijdrage level kstart+1 aan de afgeleide van th
    // for (int j=grid->jstart; j<grid->jend; ++j)
    //     #pragma ivdep
    //     for (int i=grid->istart; i<grid->iend; ++i)
    //     {
    //         const int ij  = i + j*jj;
    //         const int ijk = i + j*jj + kstart*kk;
    //         // TODO use the thermal expansion coefficient from the input later, what to do if there is no buoyancy?
    //         // Add the buoyancy production to the TKE
    //         RitPrratio = -bfluxbot[ij]/(Constants::kappa*z[kstart]*ustar[ij])*most::phih(z[kstart]/obuk[ij]) / evisc[ijk] / tPr;
    //         RitPrratio = std::min(RitPrratio, 1.-Constants::dsmall);
    //         evisc[ijk] = fac * std::sqrt(evisc[ijk]) * std::sqrt(1.-RitPrratio);
    //     }

    // NOOT: in DALES wordt K = cm*lambda*sqrt(e) ook gebruikt voor eerste gridpunt boven de bodem (gewoon gekopieerd hier)

    if(resolved_wall)
    {
        for (int k=grid->kstart; k<grid->kend; ++k)
        {
            // calculate geometric mean of filter mesh size based on Deardorff, 1973
            mlen0 = std::pow(dx*dy*dz[k], 1./3.);

            for (int j=grid->jstart; j<grid->jend; ++j)
                #pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    // Calculate eddy viscosity for momentum based on Deardorff, 1980
                    mlen = cn * std::sqrt(sgstke[ijk]) / std::sqrt(std::abs(N2[ijk]));
                    fac  = std::min(mlen0, mlen);
                    evisc[ijk] = cm * fac * std::sqrt(sgstke[ijk]) + mvisc;
                }
        }

        grid->boundary_cyclic(evisc);

        // For a resolved wall the viscosity at the wall is needed. For now, assume that the eddy viscosity
        // is zero, so set ghost cell such that the viscosity interpolated to the surface equals the molecular viscosity.
        const int kb = grid->kstart;
        const int kt = grid->kend-1;
        for (int j=0; j<grid->jcells; ++j)
            #pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                const int ijkb = i + j*jj + kb*kk;
                const int ijkt = i + j*jj + kt*kk;
                evisc[ijkb-kk] = 2 * mvisc - evisc[ijkb];
                evisc[ijkt+kk] = 2 * mvisc - evisc[ijkt];
            }
    }
    else
    {
        for (int k=grid->kstart; k<grid->kend; ++k)
        {
            // calculate geometric mean of filter mesh size based on Deardorff, 1973
            mlen0 = std::pow(dx*dy*dz[k], 1./3.);

            for (int j=grid->jstart; j<grid->jend; ++j)
                #pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    // Calculate eddy viscosity for momentum based on Deardorff, 1980
                    mlen = cn * std::sqrt(sgstke[ijk]) / std::sqrt(std::abs(N2[ijk]));
                    fac  = std::min(mlen0, mlen);
                    evisc[ijk] = cm * fac * std::sqrt(sgstke[ijk]) + mvisc;
                }
        }

        grid->boundary_cyclic(evisc);
    }
}

// NOOT : de veranderingen in deze functie zijn temptative: vrijwel identiek aan standaard calc_evisc, maar overal expliciet N2 verwijderd (vermijd delen door nul)
template <bool resolved_wall>
void Diff_sgs_tke::calc_evisc_neutral(double* restrict evisc,
                                     double* restrict sgstke, double* restrict u, double* restrict v, double* restrict w,
                                     double* restrict ufluxbot, double* restrict vfluxbot,
                                     double* restrict z, double* restrict dz, const double z0m, const double mvisc)
{
    // Variables for the wall damping.
    double mlen,mlen0,fac;

    const int jj = grid->icells;
    const int kk = grid->ijcells;

    // Make local copies to aid vectorization.
    const double dx = grid->dx;
    const double dy = grid->dy;

    const double cm  = this->cm;
    const double cn  = this->cn;

    if (resolved_wall)
    {
        for (int k=grid->kstart; k<grid->kend; ++k)
        {
            // calculate geometric mean of filter mesh size based on Deardorff, 1973
            mlen0 = std::pow(dx*dy*dz[k], 1./3.);

            for (int j=grid->jstart; j<grid->jend; ++j)
                #pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    // Calculate eddy viscosity for momentum based on Deardorff, 1980
                    mlen = cn * std::sqrt(sgstke[ijk]);
                    fac  = std::min(mlen0, mlen);
                    // NOOT : waarom wordt hier expliciet molecular viscosity bij opgeteld? En niet bijvoorbeeld bij standaard calc_evisc?
                    evisc[ijk] = cm * fac * std::sqrt(sgstke[ijk]) + mvisc;
                }
        }

        grid->boundary_cyclic(evisc);

        // For a resolved wall the viscosity at the wall is needed. For now, assume that the eddy viscosity
        // is zero, so set ghost cell such that the viscosity interpolated to the surface equals the molecular viscosity.
        const int kb = grid->kstart;
        const int kt = grid->kend-1;
        for (int j=0; j<grid->jcells; ++j)
            #pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                const int ijkb = i + j*jj + kb*kk;
                const int ijkt = i + j*jj + kt*kk;
                evisc[ijkb-kk] = 2 * mvisc - evisc[ijkb];
                evisc[ijkt+kk] = 2 * mvisc - evisc[ijkt];
            }
    }
    else
    {
        for (int k=grid->kstart; k<grid->kend; ++k)
        {
            // calculate geometric mean of filter mesh size based on Deardorff, 1973
            mlen0 = std::pow(dx*dy*dz[k], 1./3.);

            for (int j=grid->jstart; j<grid->jend; ++j)
                #pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    // Calculate eddy viscosity for momentum based on Deardorff, 1980
                    mlen = cn * std::sqrt(sgstke[ijk]);
                    fac  = std::min(mlen0, mlen);
                    // NOOT : waarom wordt hier dan NIET expliciet molecular viscosity bij opgeteld?
                    evisc[ijk] = cm * fac * std::sqrt(sgstke[ijk]);
                }
        }

        grid->boundary_cyclic(evisc);
    }

}

template <bool resolved_wall>
void Diff_sgs_tke::diff_u(double* restrict ut, double* restrict u, double* restrict v, double* restrict w,
                         double* restrict dzi, double* restrict dzhi, double* restrict evisc,
                         double* restrict fluxbot, double* restrict fluxtop,
                         double* restrict rhoref, double* restrict rhorefh)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    double eviscn, eviscs, eviscb, evisct;

    const int k_offset = resolved_wall ? 0 : 1;

    if(!resolved_wall)
    {
        // bottom boundary
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;
                eviscn = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
                eviscs = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
                evisct = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
                eviscb = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);

                ut[ijk] +=
                         // du/dx + du/dx
                         + ( evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
                           - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
                         // du/dy + dv/dx
                         + ( eviscn*((u[ijk+jj]-u[ijk   ])*dyi + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                           - eviscs*((u[ijk   ]-u[ijk-jj])*dyi + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                         // du/dz + dw/dx
                         + ( rhorefh[kstart+1] * evisct*((u[ijk+kk]-u[ijk   ])* dzhi[kstart+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
                           + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
            }

        // top boundary
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + (kend-1)*kk;
                eviscn = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
                eviscs = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
                evisct = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
                eviscb = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);
                ut[ijk] +=
                         // du/dx + du/dx
                         + ( evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
                           - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
                         // du/dy + dv/dx
                         + ( eviscn*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                           - eviscs*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                         // du/dz + dw/dx
                         + (- rhorefh[kend  ] * fluxtop[ij]
                            - rhorefh[kend-1] * eviscb*((u[ijk   ]-u[ijk-kk])* dzhi[kend-1] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) / rhoref[kend-1] * dzi[kend-1];
            }
    }

    for (int k=grid->kstart+k_offset; k<grid->kend-k_offset; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                eviscn = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
                eviscs = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
                evisct = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
                eviscb = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);
                ut[ijk] +=
                         // du/dx + du/dx
                         + ( evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
                           - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
                         // du/dy + dv/dx
                         + ( eviscn*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                           - eviscs*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                         // du/dz + dw/dx
                         + ( rhorefh[k+1] * evisct*((u[ijk+kk]-u[ijk   ])* dzhi[k+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
                           - rhorefh[k  ] * eviscb*((u[ijk   ]-u[ijk-kk])* dzhi[k  ] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) / rhoref[k] * dzi[k];
            }
}

template <bool resolved_wall>
void Diff_sgs_tke::diff_v(double* restrict vt, double* restrict u, double* restrict v, double* restrict w,
                         double* restrict dzi, double* restrict dzhi, double* restrict evisc,
                         double* restrict fluxbot, double* restrict fluxtop,
                         double* restrict rhoref, double* restrict rhorefh)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    double evisce, eviscw, eviscb, evisct;

    const int k_offset = resolved_wall ? 0 : 1;

    if(!resolved_wall)
    {
        // bottom boundary
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;
                evisce = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
                eviscw = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
                evisct = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
                eviscb = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);
                vt[ijk] +=
                         // dv/dx + du/dy
                         + ( evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                           - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                         // dv/dy + dv/dy
                         + ( evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
                           - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
                         // dv/dz + dw/dy
                         + ( rhorefh[kstart+1] * evisct*((v[ijk+kk]-v[ijk   ])*dzhi[kstart+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
                           + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
            }

        // top boundary
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + (kend-1)*kk;
                evisce = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
                eviscw = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
                evisct = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
                eviscb = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);
                vt[ijk] +=
                         // dv/dx + du/dy
                         + ( evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                           - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                         // dv/dy + dv/dy
                         + ( evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
                           - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
                         // dv/dz + dw/dy
                         + (- rhorefh[kend  ] * fluxtop[ij]
                            - rhorefh[kend-1] * eviscb*((v[ijk   ]-v[ijk-kk])*dzhi[kend-1] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) / rhoref[kend-1] * dzi[kend-1];
            }
    }

    for (int k=grid->kstart+k_offset; k<grid->kend-k_offset; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                evisce = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
                eviscw = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
                evisct = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
                eviscb = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);
                vt[ijk] +=
                         // dv/dx + du/dy
                         + ( evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                           - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                         // dv/dy + dv/dy
                         + ( evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
                           - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
                         // dv/dz + dw/dy
                         + ( rhorefh[k+1] * evisct*((v[ijk+kk]-v[ijk   ])*dzhi[k+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
                           - rhorefh[k  ] * eviscb*((v[ijk   ]-v[ijk-kk])*dzhi[k  ] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) / rhoref[k] * dzi[k];
            }
}

// NOOT :template <bool resolved_wall> // HIER OOK NODIG? Mijn 'gok' is nee, in 2de orde wordt 'afgedwongen' dat w =0 op de bodem door de term simpelweg nergens mee te nemen (?)
void Diff_sgs_tke::diff_w(double* restrict wt, double* restrict u, double* restrict v, double* restrict w,
                         double* restrict dzi, double* restrict dzhi, double* restrict evisc,
                         double* restrict rhoref, double* restrict rhorefh)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    double evisce, eviscw, eviscn, eviscs;

    for (int k=grid->kstart+1; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                evisce = 0.25*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+ii-kk] + evisc[ijk+ii]);
                eviscw = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-ii] + evisc[ijk   -kk] + evisc[ijk   ]);
                eviscn = 0.25*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+jj-kk] + evisc[ijk+jj]);
                eviscs = 0.25*(evisc[ijk-jj-kk] + evisc[ijk-jj] + evisc[ijk   -kk] + evisc[ijk   ]);
                wt[ijk] +=
                         // dw/dx + du/dz
                         + ( evisce*((w[ijk+ii]-w[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-kk])*dzhi[k])
                           - eviscw*((w[ijk   ]-w[ijk-ii])*dxi + (u[ijk   ]-u[ijk+  -kk])*dzhi[k]) ) * dxi
                         // dw/dy + dv/dz
                         + ( eviscn*((w[ijk+jj]-w[ijk   ])*dyi + (v[ijk+jj]-v[ijk+jj-kk])*dzhi[k])
                           - eviscs*((w[ijk   ]-w[ijk-jj])*dyi + (v[ijk   ]-v[ijk+  -kk])*dzhi[k]) ) * dyi
                         // dw/dz + dw/dz
                         + ( rhoref[k  ] * evisc[ijk   ]*(w[ijk+kk]-w[ijk   ])*dzi[k  ]
                           - rhoref[k-1] * evisc[ijk-kk]*(w[ijk   ]-w[ijk-kk])*dzi[k-1] ) / rhorefh[k] * 2.* dzhi[k];
            }
}

template <bool resolved_wall>
void Diff_sgs_tke::diff_c(double* restrict at, double* restrict a,
                         double* restrict dz, double* restrict dzi, double* restrict dzhi, double* restrict evisc, double* restrict sgstke, double* restrict N2,
                         double* restrict fluxbot, double* restrict fluxtop,
                         double* restrict rhoref, double* restrict rhorefh)
{
    // Variables for the length scale
    double mlen,mlen0,fac;
    // Variables for the stability dependent turbulent Prandtl number
    double ch,tPri;

    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    const double dx     = grid->dx;
    const double dy     = grid->dy;
    const double dxidxi = 1./(grid->dx * grid->dx);
    const double dyidyi = 1./(grid->dy * grid->dy);

    // local copies to aid vectorization
    // NOOT: volgens mij is het maken van local copies onnodig voor consts ?
    double ch1 = this->ch1;
    double ch2 = this->ch2;
    double cn  = this->cn;

    double evisce,eviscw,eviscn,eviscs,evisct,eviscb;

    const int k_offset = resolved_wall ? 0 : 1;

    if(!resolved_wall)
    {
        // bottom boundary
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;

                // calculate geometric mean of filter mesh size based on Deardorff, 1973
                mlen0 = std::pow(dx*dy*dz[kstart], 1./3.);

                // Calculate turbulent length scale based on Deardorff, 1980
                mlen = cn * std::sqrt(sgstke[ijk]) / std::sqrt(std::abs(N2[ijk]));
                fac  = std::min(mlen0, mlen);

                // Calculate the inverse stability dependent turbulent Prandtl number
                tPri = (ch1 + ch2 * fac / mlen0);

                evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii]) * tPri;
                eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ]) * tPri;
                eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj]) * tPri;
                eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ]) * tPri;
                evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk]) * tPri;
                eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ]) * tPri;

                at[ijk] +=
                         + ( evisce*(a[ijk+ii]-a[ijk   ])
                           - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi
                         + ( eviscn*(a[ijk+jj]-a[ijk   ])
                           - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                         + ( rhorefh[kstart+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[kstart+1]
                           + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
            }

        // top boundary
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + (kend-1)*kk;

                // calculate geometric mean of filter mesh size based on Deardorff, 1973
                mlen0 = std::pow(dx*dy*dz[kend-1], 1./3.);

                // Calculate turbulent length scale based on Deardorff, 1980
                mlen = cn * std::sqrt(sgstke[ijk]) / std::sqrt(std::abs(N2[ijk]));
                fac  = std::min(mlen0, mlen);

                // Calculate the inverse stability dependent turbulent Prandtl number
                tPri = (ch1 + ch2 * fac / mlen0);

                evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii]) * tPri;
                eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ]) * tPri;
                eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj]) * tPri;
                eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ]) * tPri;
                evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk]) * tPri;
                eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ]) * tPri;

                at[ijk] +=
                         + ( evisce*(a[ijk+ii]-a[ijk   ])
                           - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi
                         + ( eviscn*(a[ijk+jj]-a[ijk   ])
                           - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                         + (-rhorefh[kend  ] * fluxtop[ij]
                           - rhorefh[kend-1] * eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[kend-1] ) / rhoref[kend-1] * dzi[kend-1];
            }
    }

    for (int k=grid->kstart+k_offset; k<grid->kend-k_offset; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // calculate geometric mean of filter mesh size based on Deardorff, 1973
                mlen0 = std::pow(dx*dy*dz[k], 1./3.);

                // Calculate turbulent length scale based on Deardorff, 1980
                mlen = cn * std::sqrt(sgstke[ijk]) / std::sqrt(std::abs(N2[ijk]));
                fac  = std::min(mlen0, mlen);

                // Calculate the inverse stability dependent turbulent Prandtl number
                tPri = (ch1 + ch2 * fac / mlen0);

                evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii]) * tPri;
                eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ]) * tPri;
                eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj]) * tPri;
                eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ]) * tPri;
                evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk]) * tPri;
                eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ]) * tPri;

                at[ijk] +=
                         + ( evisce*(a[ijk+ii]-a[ijk   ])
                           - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi
                         + ( eviscn*(a[ijk+jj]-a[ijk   ])
                           - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                         + ( rhorefh[k+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[k+1]
                           - rhorefh[k  ] * eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[k]  ) / rhoref[k] * dzi[k];
            }
}

template <bool resolved_wall> // Maak de resolved wall optie nog
void Diff_sgs_tke::diff_sgs_tke(double* restrict sgstket, double* restrict sgstke,
                         double* restrict dzi, double* restrict dzhi, double* restrict evisc,
                         double* restrict fluxbot, double* restrict fluxtop,
                         double* restrict rhoref, double* restrict rhorefh)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    const double dxidxi = 1./(grid->dx * grid->dx);
    const double dyidyi = 1./(grid->dy * grid->dy);

    double evisce,eviscw,eviscn,eviscs,evisct,eviscb;

    const int k_offset = resolved_wall ? 0 : 1;

    if(!resolved_wall)
    {
        // bottom boundary
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;
                evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii]);
                eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ]);
                eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj]);
                eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ]);
                evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk]);
                eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ]);

                sgstket[ijk] +=
                             + ( evisce*(sgstke[ijk+ii]-sgstke[ijk   ])
                               - eviscw*(sgstke[ijk   ]-sgstke[ijk-ii]) ) * dxidxi
                             + ( eviscn*(sgstke[ijk+jj]-sgstke[ijk   ])
                               - eviscs*(sgstke[ijk   ]-sgstke[ijk-jj]) ) * dyidyi
                             + ( rhorefh[kstart+1] * evisct*(sgstke[ijk+kk]-sgstke[ijk   ])*dzhi[kstart+1]
                               ) / rhoref[kstart] * dzi[kstart];
                               // NOOT : Voor nu ervoor gekozen om diffusie via fluxbot weg te halen
                            // + ( rhorefh[kstart+1] * evisct*(sgstke[ijk+kk]-sgstke[ijk   ])*dzhi[kstart+1]
                            //   + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
            }

        // top boundary
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + (kend-1)*kk;
                evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii]);
                eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ]);
                eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj]);
                eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ]);
                evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk]);
                eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ]);

                sgstket[ijk] +=
                             + ( evisce*(sgstke[ijk+ii]-sgstke[ijk   ])
                               - eviscw*(sgstke[ijk   ]-sgstke[ijk-ii]) ) * dxidxi
                             + ( eviscn*(sgstke[ijk+jj]-sgstke[ijk   ])
                               - eviscs*(sgstke[ijk   ]-sgstke[ijk-jj]) ) * dyidyi
                             + ( - rhorefh[kend-1] * eviscb*(sgstke[ijk   ]-sgstke[ijk-kk])*dzhi[kend-1] ) / rhoref[kend-1] * dzi[kend-1];
                             // NOOT: voor nu ervoor gekozen om diffusie door top via fluxtop weg te halen
                             //+ (-rhorefh[kend  ] * fluxtop[ij]
                            //   - rhorefh[kend-1] * eviscb*(sgstke[ijk   ]-sgstke[ijk-kk])*dzhi[kend-1] ) / rhoref[kend-1] * dzi[kend-1];
            }

    }

    // NOOT : Resolved wall voor SGS TKE vereist dat de ghostcell wordt ingesteld -> geen flux vereist dat (sgstke[ijk   ]-sgstke[ijk-kk]) == 0 --> dus afgeleide is nul --> swbot[sgs_tke] = neumann, sbot[sgs_tke] = 0
    // NOOT : voor 'veiligheid' kan dit hardcoded worden, of een error/warning geprint worden indien het niet correct ingesteld is ??
    for (int k=grid->kstart+k_offset; k<grid->kend-k_offset; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii]);
                eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ]);
                eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj]);
                eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ]);
                evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk]);
                eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ]);

                sgstket[ijk] +=
                             + ( evisce*(sgstke[ijk+ii]-sgstke[ijk   ])
                               - eviscw*(sgstke[ijk   ]-sgstke[ijk-ii]) ) * dxidxi
                             + ( eviscn*(sgstke[ijk+jj]-sgstke[ijk   ])
                               - eviscs*(sgstke[ijk   ]-sgstke[ijk-jj]) ) * dyidyi
                             + ( rhorefh[k+1] * evisct*(sgstke[ijk+kk]-sgstke[ijk   ])*dzhi[k+1]
                               - rhorefh[k  ] * eviscb*(sgstke[ijk   ]-sgstke[ijk-kk])*dzhi[k]  ) / rhoref[k] * dzi[k];
            }

}

void Diff_sgs_tke::calc_sgs_tke_shear_tend_2nd(double* restrict sgstket,
                                              double* restrict evisc, double* restrict strain2)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    // NOOT: strain2 en eddy viscosity zijn gedefinieerd op cell center niveau? Dan zou productie direct moeten volgen uit hun product
    for (int k=grid->kstart; k<grid->kend; ++k)
    {
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                // Calculate shear production of TKE based on Deardorff, 1980
                sgstket[ijk] = evisc[ijk] * strain2[ijk];
            }
    }
}

void Diff_sgs_tke::calc_sgs_tke_buoyancy_tend_2nd(double* restrict sgstket, double* restrict sgstke,
                                                 double* restrict evisc, double* restrict N2, double* restrict dz)
{
    // NOOT: N2 en eddy viscosity zijn gedefinieerd op cell center niveau? Dan zou productie direct moeten volgen uit hun product

    // Variables for the length scale
    double mlen,mlen0,fac;
    // Variables for the stability dependent turbulent Prandtl number
    double ch,tPri;

    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const double dx = grid->dx;
    const double dy = grid->dy;

    // local copies to aid vectorization
    // NOOT: volgens mij is het maken van local copies onnodig voor consts ?
    double ch1 = this->ch1;
    double ch2 = this->ch2;
    double cn  = this->cn;

    for (int k=grid->kstart; k<grid->kend; ++k)
    {
        // calculate geometric mean of filter mesh size based on Deardorff, 1973
        mlen0 = std::pow(dx*dy*dz[k], 1./3.);

        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // Calculate turbulent length scale based on Deardorff, 1980
                mlen = cn * std::sqrt(sgstke[ijk]) / std::sqrt(std::abs(N2[ijk]));
                fac  = std::min(mlen0, mlen);

                // Calculate the inverse stability dependent turbulent Prandtl number
                tPri = (ch1 + ch2 * fac / mlen0);

                // Calculate buoyancy production of TKE based on Deardorff, 1980
                sgstket[ijk] = - evisc[ijk] * N2[ijk] * tPri;
            }
    }
}

// NOOT: voor dissipatie is geen 2de en 4de orde nodig --> lokaal, in grid cell effect
// NOOT: N2 wordt nu berekend via 2de orde interpolatie --> dus zelfs wanneer 4de orde schema voor LES ingeprogrammeerd wordt, blijft geheel toch lagere orde?
void Diff_sgs_tke::calc_sgs_tke_dissipation_2nd(double* restrict sgstket, double* restrict sgstke,
                                               double* restrict N2, double* restrict dz)
{
    // Variables for the length scale and dissipation.
    double mlen,mlen0,fac;
    double ce;

    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    const double dx = grid->dx;
    const double dy = grid->dy;

    // Local copies to aid vectorization
    double ce1 = this->ce1;
    double ce2 = this->ce2;
    double cn  = this->cn;

    for (int k=grid->kstart; k<grid->kend; ++k)
    {
        // calculate geometric mean of filter mesh size based on Deardorff, 1973
        mlen0 = std::pow(dx*dy*dz[k], 1./3.);

        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                // Calculate dissipation of TKE based on Deardorff, 1980
                mlen         = cn * std::sqrt(sgstke[ijk]) / std::sqrt(std::abs(N2[ijk]));
                fac          = std::min(mlen0, mlen);
                ce           = ce1 + ce2 * fac / mlen0;
                sgstket[ijk] = -1 * ce * std::pow(sgstke[ijk], 3./2.) / fac;
            }
    }
}

//void Diff_sgs_tke::calc_sgs_tke_dissipation_4th(double* restrict at, double* restrict a,)
//{}
//void Diff_sgs_tke::calc_sgs_tke_buoyancy_tend_4th(double* restrict at, double* restrict a,)
//{}
//void Diff_sgs_tke::calc_sgs_tke_shear_tend_4th(double* restrict at, double* restrict a,)
//{}

// NOOT: turbulent Prandtl getal hier moet dus variabel worden, i.e., hoogte- en stabiliteitsafhankelijk
// om te vermijden het gehele turbulente Prandtl veld op te slaan, hier meegeven tPr = 2 ?
// dit kan de code in theorie trager maken wanneer dit de beperkende tijdstap wordt. Voor SBLs eigenlijk geen probleem
double Diff_sgs_tke::calc_dnmul(double* restrict evisc, double* restrict dzi, double tPr)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const double dxidxi = 1./(grid->dx * grid->dx);
    const double dyidyi = 1./(grid->dy * grid->dy);

    // NOOT: vanwaar dit statement? Het is mogelijk dat tPr groter wordt dan 1. In dat geval zou dit statement ongewenst zijn?
    const double tPrfac = std::min(1., tPr);
    double dnmul = 0;

    // get the maximum time step for diffusion
    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                dnmul = std::max(dnmul, std::abs(tPrfac*evisc[ijk]*(dxidxi + dyidyi + dzi[k]*dzi[k])));
            }

    grid->get_max(&dnmul);

    return dnmul;
}

void Diff_sgs_tke::calc_prandtl(double* restrict prandtl, double* restrict sgstke, double* restrict N2, double* restrict dz)
{

      // Variables for the length scale
      double mlen,mlen0,fac;
      // Variables for the stability dependent turbulent Prandtl number
      double ch;

      const int jj = grid->icells;
      const int kk = grid->ijcells;

      const double dx = grid->dx;
      const double dy = grid->dy;

      // local copies to aid vectorization
      // NOOT: volgens mij is het maken van local copies onnodig voor consts ?
      double ch1 = this->ch1;
      double ch2 = this->ch2;
      double cn  = this->cn;

      for (int k=grid->kstart; k<grid->kend; ++k)
      {
          // calculate geometric mean of filter mesh size based on Deardorff, 1973
          mlen0 = std::pow(dx*dy*dz[k], 1./3.);

          for (int j=grid->jstart; j<grid->jend; ++j)
              #pragma ivdep
              for (int i=grid->istart; i<grid->iend; ++i)
              {
                  const int ijk = i + j*jj + k*kk;

                  // Calculate turbulent length scale based on Deardorff, 1980
                  mlen = cn * std::sqrt(sgstke[ijk]) / std::sqrt(std::abs(N2[ijk]));
                  fac  = std::min(mlen0, mlen);

                  // Calculate the stability dependent turbulent Prandtl number
                  prandtl[ijk] = 1 / (ch1 + ch2 * fac / mlen0);
              }
      }
}

void Diff_sgs_tke::init_stats() ///< function for additional statistics, SvdLinden, May 2018
{
    if (model->stats->get_switch() == "1")
    {
        model->stats->add_prof("tPr", "Turbulent Prandtl number", "-", "z"); // stability dependent turbulent Prandtl number
    }
}

void Diff_sgs_tke::exec_stats(Mask *m)
{

    if (model->thermo->get_switch() != "0") // actually only do statistics of Prandtl number if there is one
    {
        const double NoOffset = 0.;

        // define the location
        const int sloc[] = {0,0,0};

        // retrieve buoyancy field for use in SGS TKE production and store in tmp1
        // store the buoyancyflux in tmp1
        model->thermo->get_buoyancy_fluxbot(fields->atmp["tmp1"]);
        // retrieve the full field in tmp1 and use tmp2 for temporary calculations
        model->thermo->get_thermo_field(fields->atmp["tmp1"], fields->atmp["tmp2"], "N2", false);

        // recalculate stability dependent turbulent Prandtl number and store in tmp2
        calc_prandtl(fields->atmp["tmp2"]->data, fields->sp["sgstke"]->data, fields->atmp["tmp1"]->data, grid->dz);

        // calculate the mean, standard mask is always stored in tmp3 (see model.cxx)
        model->stats->calc_mean(m->profs["tPr"].data, fields->atmp["tmp2"]->data, NoOffset, sloc,
              fields->atmp["tmp3"]->data, model->stats->nmask);
    }
}
