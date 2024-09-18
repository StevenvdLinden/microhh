/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
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

#include "immersed_boundary.h"
#include "boundary.h"
#include "fields.h"
#include "tools.h"
#include "fast_math.h"

namespace
{
    namespace fm = Fast_math;
    
    template<typename TF> __global__
    void set_forcing_points_u_g(
        TF* const __restrict__ tend_u,
        const TF* const __restrict__ tend_v,
        const TF* const __restrict__ tend_w,
        const TF* const __restrict__ fld_u,
        const TF* const __restrict__ fld_v,
        const TF* const __restrict__ fld_w,
        const TF* const __restrict__ boundary_value,
        const int* const __restrict__ gi, const int* const __restrict__ gj, const int* const __restrict__ gk,
        const TF* const __restrict__ rot,
        const int* const __restrict__ ipui, const int* const __restrict__ ipuj, const int* const __restrict__ ipuk, const TF* const __restrict__ c_idw_u,
        const int* const __restrict__ ipvi, const int* const __restrict__ ipvj, const int* const __restrict__ ipvk, const TF* const __restrict__ c_idw_v, 
        const int* const __restrict__ ipwi, const int* const __restrict__ ipwj, const int* const __restrict__ ipwk, const TF* const __restrict__ c_idw_w,
        const int* const __restrict__ ipsi, const int* const __restrict__ ipsj, const int* const __restrict__ ipsk, const TF* const __restrict__ c_idw_s, // SvdL, 20240901: not used for now..
        const TF* const __restrict__ db, const TF* const __restrict__ di, const TF* const __restrict__ z0b,
        Boundary_type bc, const TF visc, const int n_fpoints, const int n_idw_loc,
        const int icells, const int ijcells,
        const double dt)
    {
        const int n    = blockIdx.x*blockDim.x + threadIdx.x;

        const int rdim = 9;                                        
        const TF  dtf  = static_cast<TF>(dt); // SvdL, 20240901: seems like ugly solution, see how to improve later.
 
        TF u_ip_la;
        TF v_ip_la;
        TF w_ip_la;
        TF u_fp_la;
        TF v_fp_la;
        TF w_fp_la;

        if (n < n_fpoints)
        {
            const int ijkf = gi[n] + gj[n] * icells + gk[n] * ijcells; // field location of forcing point
            const TF r11 = rot[rdim * n];                              // this is maybe a redundant type defintion, since incoming rot is already type <TF>
            const TF r12 = rot[rdim * n + 1];
            const TF r13 = rot[rdim * n + 2];
            const TF r21 = rot[rdim * n + 3];
            const TF r22 = rot[rdim * n + 4];
            const TF r23 = rot[rdim * n + 5];
            const TF r31 = rot[rdim * n + 6];
            const TF r32 = rot[rdim * n + 7];
            const TF r33 = rot[rdim * n + 8];

            TF u_ip = TF(0.);
            TF v_ip = TF(0.);
            TF w_ip = TF(0.);

            // 1. interpolate surroundings neighbours to interpolation point
            for (int i = 0; i < n_idw_loc; ++i)
            {
                const int ii = i + n * n_idw_loc;
                const int ijku = ipui[ii] + ipuj[ii] * icells + ipuk[ii] * ijcells;
                const int ijkv = ipvi[ii] + ipvj[ii] * icells + ipvk[ii] * ijcells;
                const int ijkw = ipwi[ii] + ipwj[ii] * icells + ipwk[ii] * ijcells;
                
                // Do the correction based on the auxiliary velocity (i.e. intermediate velocity at next timestep without pressure forcing).
                u_ip += c_idw_u[ii] * (fld_u[ijku] + dtf * tend_u[ijku] );
                v_ip += c_idw_v[ii] * (fld_v[ijkv] + dtf * tend_v[ijkv] );
                w_ip += c_idw_w[ii] * (fld_w[ijkw] + dtf * tend_w[ijkw] );
            }

            // 2. rotate velocities to locally align with surface tangent (under the assumption that flow at second layer still aligns)
            u_ip_la = r11 * u_ip + r12 * v_ip + r13 * w_ip;
            v_ip_la = r21 * u_ip + r22 * v_ip + r23 * w_ip;
            w_ip_la = r31 * u_ip + r32 * v_ip + r33 * w_ip;

            // for now, (1) neglect flow rotation over height, (2) neglect stability effects (requires "fine enough" grid),
            // (3) assume both points are in logarithmic layer, and (4) assume zero-valued Dirichlet conditions for momentum (i.e. no-slip condition)
            // future options: investigate use of Van Driest (1956) correction and/or DNS mode.
            if (db[n] > z0b[n])
            {
                // 3. calculate (locally-aligned) velocity at forcing point
                u_fp_la = u_ip_la * std::log(db[n] / z0b[n]) / std::log(di[n] / z0b[n]);
                v_fp_la = v_ip_la * std::log(db[n] / z0b[n]) / std::log(di[n] / z0b[n]);
                w_fp_la = w_ip_la * fm::pow2(db[n] / di[n]);

                // 4. rotate back to standard grid alginment (only one component is needed here), 
                // AND overwrite old tendency at forcing point with new one to achieve this.
                tend_u[ijkf] = ( (r11 * u_fp_la + r21 * v_fp_la + r31 * w_ip_la) - fld_u[ijkf] ) / dtf;
            }
            else // SvdL, 20240901: investigate if can be changed into Van Driest like correction..
            {
                tend_u[ijkf] = ( TF(0.) - fld_u[ijkf] ) / dtf;
            }

        }
    }

    template<typename TF> __global__
    void set_forcing_points_v_g(
        const TF* const __restrict__ tend_u,
        TF* const __restrict__ tend_v,
        const TF* const __restrict__ tend_w,
        const TF* const __restrict__ fld_u,
        const TF* const __restrict__ fld_v,
        const TF* const __restrict__ fld_w,
        const TF* const __restrict__ boundary_value,
        const int* const __restrict__ gi, const int* const __restrict__ gj, const int* const __restrict__ gk,
        const TF* const __restrict__ rot,
        const int* const __restrict__ ipui, const int* const __restrict__ ipuj, const int* const __restrict__ ipuk, const TF* const __restrict__ c_idw_u,
        const int* const __restrict__ ipvi, const int* const __restrict__ ipvj, const int* const __restrict__ ipvk, const TF* const __restrict__ c_idw_v, 
        const int* const __restrict__ ipwi, const int* const __restrict__ ipwj, const int* const __restrict__ ipwk, const TF* const __restrict__ c_idw_w,
        const int* const __restrict__ ipsi, const int* const __restrict__ ipsj, const int* const __restrict__ ipsk, const TF* const __restrict__ c_idw_s, 
        const TF* const __restrict__ db, const TF* const __restrict__ di, const TF* const __restrict__ z0b,
        Boundary_type bc, const TF visc, const int n_fpoints, const int n_idw_loc,
        const int icells, const int ijcells,
        const double dt)
    {
        const int n    = blockIdx.x*blockDim.x + threadIdx.x;

        const int rdim = 9;                                       
        const TF  dtf  = static_cast<TF>(dt); // SvdL, 20240901: seems like ugly solution, see how to improve later.

        TF u_ip_la;
        TF v_ip_la;
        TF w_ip_la;
        TF u_fp_la;
        TF v_fp_la;
        TF w_fp_la;

        // Loop over all points to be forced
        if (n < n_fpoints)
        {
            const int ijkf = gi[n] + gj[n] * icells + gk[n] * ijcells; // field location of forcing point
            const TF r11 = rot[rdim * n];
            const TF r12 = rot[rdim * n + 1];
            const TF r13 = rot[rdim * n + 2];
            const TF r21 = rot[rdim * n + 3];
            const TF r22 = rot[rdim * n + 4];
            const TF r23 = rot[rdim * n + 5];
            const TF r31 = rot[rdim * n + 6];
            const TF r32 = rot[rdim * n + 7];
            const TF r33 = rot[rdim * n + 8];

            TF u_ip = TF(0.);
            TF v_ip = TF(0.);
            TF w_ip = TF(0.);

            // 1. interpolate surroundings neighbours to interpolation point
            for (int i = 0; i < n_idw_loc; ++i)
            {
                const int ii = i + n * n_idw_loc;                                  
                const int ijku = ipui[ii] + ipuj[ii] * icells + ipuk[ii] * ijcells;
                const int ijkv = ipvi[ii] + ipvj[ii] * icells + ipvk[ii] * ijcells;
                const int ijkw = ipwi[ii] + ipwj[ii] * icells + ipwk[ii] * ijcells;

                // Do the correction based on the auxiliary velocity (i.e. intermediate velocity at next timestep without pressure forcing).
                u_ip += c_idw_u[ii] * (fld_u[ijku] + dtf * tend_u[ijku] );
                v_ip += c_idw_v[ii] * (fld_v[ijkv] + dtf * tend_v[ijkv] );
                w_ip += c_idw_w[ii] * (fld_w[ijkw] + dtf * tend_w[ijkw] );
            }

            // 2. rotate velocities to locally align with surface tangent (under the assumption that flow at second layer still aligns)
            u_ip_la = r11 * u_ip + r12 * v_ip + r13 * w_ip;
            v_ip_la = r21 * u_ip + r22 * v_ip + r23 * w_ip;
            w_ip_la = r31 * u_ip + r32 * v_ip + r33 * w_ip;

            // for now, (1) neglect flow rotation over height, (2) neglect stability effects (requires "fine enough" grid),
            // (3) assume both points are in logarithmic layer, and (4) assume zero-valued Dirichlet conditions for momentum (i.e. no-slip condition)
            // future options: investigate use of Van Driest (1956) correction and/or DNS mode.
            if (db[n] > z0b[n])
            {
                // 3. calculate (locally-aligned) velocity at forcing point
                u_fp_la = u_ip_la * std::log(db[n] / z0b[n]) / std::log(di[n] / z0b[n]);
                v_fp_la = v_ip_la * std::log(db[n] / z0b[n]) / std::log(di[n] / z0b[n]);
                w_fp_la = w_ip_la *  fm::pow2(db[n] / di[n]);

                // 4. rotate back to standard grid alginment (only one component is needed here), 
                // AND overwrite old tendency at forcing point with new one to achieve this.
                tend_v[ijkf] = ( (r12 * u_fp_la + r22 * v_fp_la + r32 * w_fp_la) - fld_v[ijkf] ) / dtf;
            }
            else // SvdL, 29-06-2023: change later into Van Driest like correction..
            {
                tend_v[ijkf] = ( TF(0.) - fld_v[ijkf] ) / dtf;
            }
        }
    }

    template<typename TF> __global__
    void set_forcing_points_w_g(
        const TF* const __restrict__ tend_u,
        const TF* const __restrict__ tend_v,
        TF* const __restrict__ tend_w,
        const TF* const __restrict__ fld_u,
        const TF* const __restrict__ fld_v,
        const TF* const __restrict__ fld_w,
        const TF* const __restrict__ boundary_value,
        const int* const __restrict__ gi, const int* const __restrict__ gj, const int* const __restrict__ gk,
        const TF* const __restrict__ rot,
        const int* const __restrict__ ipui, const int* const __restrict__ ipuj, const int* const __restrict__ ipuk, const TF* const __restrict__ c_idw_u,
        const int* const __restrict__ ipvi, const int* const __restrict__ ipvj, const int* const __restrict__ ipvk, const TF* const __restrict__ c_idw_v, 
        const int* const __restrict__ ipwi, const int* const __restrict__ ipwj, const int* const __restrict__ ipwk, const TF* const __restrict__ c_idw_w,
        const int* const __restrict__ ipsi, const int* const __restrict__ ipsj, const int* const __restrict__ ipsk, const TF* const __restrict__ c_idw_s, // SvdL, 20240901: not used for now..
        const TF* const __restrict__ db, const TF* const __restrict__ di, const TF* const __restrict__ z0b,
        Boundary_type bc, const TF visc, const int n_fpoints, const int n_idw_loc,
        const int icells, const int ijcells,
        const double dt)
    {   
        const int n    = blockIdx.x*blockDim.x + threadIdx.x;

        const int rdim = 9;                                        
        const TF  dtf  = static_cast<TF>(dt); // SvdL, 20240901: seems like ugly solution, see how to improve later.

        TF u_ip_la;
        TF v_ip_la;
        TF w_ip_la;
        TF u_fp_la;
        TF v_fp_la;
        TF w_fp_la;

        // Loop over all points to be forced
        if (n < n_fpoints)
        {
            const int ijkf = gi[n] + gj[n] * icells + gk[n] * ijcells; // field location of forcing point
            const TF r11 = rot[rdim * n];
            const TF r12 = rot[rdim * n + 1];
            const TF r13 = rot[rdim * n + 2];
            const TF r21 = rot[rdim * n + 3];
            const TF r22 = rot[rdim * n + 4];
            const TF r23 = rot[rdim * n + 5];
            const TF r31 = rot[rdim * n + 6];
            const TF r32 = rot[rdim * n + 7];
            const TF r33 = rot[rdim * n + 8];

            TF u_ip = TF(0.);
            TF v_ip = TF(0.);
            TF w_ip = TF(0.);

            // 1. interpolate surroundings neighbours to interpolation point
            for (int i = 0; i < n_idw_loc; ++i)
            {
                const int ii = i + n * n_idw_loc;                                  
                const int ijku = ipui[ii] + ipuj[ii] * icells + ipuk[ii] * ijcells;
                const int ijkv = ipvi[ii] + ipvj[ii] * icells + ipvk[ii] * ijcells;
                const int ijkw = ipwi[ii] + ipwj[ii] * icells + ipwk[ii] * ijcells;

                // Do the correction based on the auxiliary velocity (i.e. intermediate velocity at next timestep without pressure forcing).
                u_ip += c_idw_u[ii] * (fld_u[ijku] + dtf * tend_u[ijku] );
                v_ip += c_idw_v[ii] * (fld_v[ijkv] + dtf * tend_v[ijkv] );
                w_ip += c_idw_w[ii] * (fld_w[ijkw] + dtf * tend_w[ijkw] );
            }

            // 2. rotate velocities to locally align with surface tangent (under the assumption that flow at second layer still aligns)
            u_ip_la = r11 * u_ip + r12 * v_ip + r13 * w_ip;
            v_ip_la = r21 * u_ip + r22 * v_ip + r23 * w_ip;
            w_ip_la = r31 * u_ip + r32 * v_ip + r33 * w_ip;

            // for now, (1) neglect flow rotation over height, (2) neglect stability effects (requires "fine enough" grid),
            // (3) assume both points are in logarithmic layer, and (4) assume zero-valued Dirichlet conditions for momentum (i.e. no-slip condition)
            // future options: investigate use of Van Driest (1956) correction and/or DNS mode.
            if (db[n] > z0b[n])
            {
                // 3. calculate (locally-aligned) velocity at forcing point
                u_fp_la = u_ip_la * std::log(db[n] / z0b[n]) / std::log(di[n] / z0b[n]);
                v_fp_la = v_ip_la * std::log(db[n] / z0b[n]) / std::log(di[n] / z0b[n]);
                w_fp_la = w_ip_la * fm::pow2(db[n] / di[n]);

                // 4. rotate back to standard grid alginment (only one component is needed here), 
                // AND overwrite old tendency at forcing point with new one to achieve this.
                tend_w[ijkf] = ( (r13 * u_fp_la + r23 * v_fp_la + r33 * w_fp_la) - fld_w[ijkf] ) / dtf;
            }
            else // SvdL, 29-06-2023: change later into Van Driest like correction..
            {
                tend_w[ijkf] = ( TF(0.) - fld_w[ijkf] ) /dtf;
            }
        }
    }

    template<typename TF> __global__
    void set_forcing_points_c_g(
        TF* const __restrict__ tend_c,
        const TF* const __restrict__ fld_c,
        const TF* const __restrict__ fld_u,
        const TF* const __restrict__ fld_v,
        const TF* const __restrict__ fld_w,
        const TF* const __restrict__ boundary_value,
        const int* const __restrict__ gi, const int* const __restrict__ gj, const int* const __restrict__ gk,
        const TF* const __restrict__ rot,
        const int* const __restrict__ ipui, const int* const __restrict__ ipuj, const int* const __restrict__ ipuk, const TF* const __restrict__ c_idw_u,
        const int* const __restrict__ ipvi, const int* const __restrict__ ipvj, const int* const __restrict__ ipvk, const TF* const __restrict__ c_idw_v, 
        const int* const __restrict__ ipwi, const int* const __restrict__ ipwj, const int* const __restrict__ ipwk, const TF* const __restrict__ c_idw_w,
        const int* const __restrict__ ipsi, const int* const __restrict__ ipsj, const int* const __restrict__ ipsk, const TF* const __restrict__ c_idw_s,
        const TF* const __restrict__ db, const TF* const __restrict__ di, const TF* const __restrict__ z0b,
        Boundary_type bc, const TF visc, const int n_fpoints, const int n_idw_loc,
        const int icells, const int ijcells, 
        const double dt)
    {
        const int n    = blockIdx.x*blockDim.x + threadIdx.x;

        const TF  dtf  = static_cast<TF>(dt); // SvdL, 20240901: seems like ugly solution, see how to improve later.

        // For Dirichlet BCs
        if (bc == Boundary_type::Dirichlet_type)
        {
            // Loop over all points to be forced
            if (n < n_fpoints)
            {
                const int ijkf = gi[n] + gj[n] * icells + gk[n] * ijcells; // field location of forcing point

                TF c_ip = TF(0.);

                // 1. interpolate surroundings neighbours to interpolation point
                for (int i = 0; i < n_idw_loc; ++i)
                {
                    const int ii = i + n * n_idw_loc;                                  
                    const int ijki = ipsi[ii] + ipsj[ii] * icells + ipsk[ii] * ijcells;
                    c_ip += c_idw_s[ii] * (fld_c[ijki] + dtf * tend_c[ijki] );
                }

                // 2. calculate scalar at forcing point and force at once
                // for now, (1) neglect stability effects (requires "fine enough" grid), and (2) assume both points are in logarithmic layer
                if (db[n] > z0b[n])
                {
                    tend_c[ijkf] = ( ( (c_ip - boundary_value[n]) * std::log(db[n] / z0b[n]) / std::log( di[n] / z0b[n]) + boundary_value[n] ) - fld_c[ijkf] ) / dtf;
                }
                else
                {
                    tend_c[ijkf] = ( boundary_value[n] - fld_c[ijkf] ) / dtf;
                }
            }
        }
        else if (bc == Boundary_type::Flux_type)
        {
            return; // SvdL, 20240918: still implement later
        }
        else if (bc == Boundary_type::Neumann_type)
        {
            return; // SvdL, 20240918: still implement later
        }
    }

    template<typename TF> __global__
    void set_forcing_points_evisc_g(
        TF* const __restrict__ fld_evisc,
        const TF* const __restrict__ fld_u,
        const TF* const __restrict__ fld_v,
        const TF* const __restrict__ fld_w,
        // const TF* const __restrict__ boundary_value, //<< SvdL, 20240909: for now not needed for eddy viscosity 
        const int* const __restrict__ gi, const int* const __restrict__ gj, const int* const __restrict__ gk,
        const TF* const __restrict__ rot,
        const int* const __restrict__ ipui, const int* const __restrict__ ipuj, const int* const __restrict__ ipuk, const TF* const __restrict__ c_idw_u,
        const int* const __restrict__ ipvi, const int* const __restrict__ ipvj, const int* const __restrict__ ipvk, const TF* const __restrict__ c_idw_v, 
        const int* const __restrict__ ipwi, const int* const __restrict__ ipwj, const int* const __restrict__ ipwk, const TF* const __restrict__ c_idw_w,
        const int* const __restrict__ ipsi, const int* const __restrict__ ipsj, const int* const __restrict__ ipsk, const TF* const __restrict__ c_idw_s,
        const TF* const __restrict__ db, const TF* const __restrict__ di, const TF* const __restrict__ z0b,
        Boundary_type bc, const TF visc, const int n_fpoints, const int n_idw_loc,
        const int icells, const int ijcells)
    {
        const int n    = blockIdx.x*blockDim.x + threadIdx.x;

        const int rdim = 9;

        TF u_ip_la;
        TF v_ip_la;
        TF w_ip_la;
        TF umag_ip_la;
        // TF ustar_ip;

        // Loop over all points to be forced
        if (n < n_fpoints)
        {
            const int ijkf = gi[n] + gj[n] * icells + gk[n] * ijcells; // field location of forcing point
            const TF r11 = rot[rdim * n];
            const TF r12 = rot[rdim * n + 1];
            const TF r13 = rot[rdim * n + 2];
            const TF r21 = rot[rdim * n + 3];
            const TF r22 = rot[rdim * n + 4];
            const TF r23 = rot[rdim * n + 5];
            const TF r31 = rot[rdim * n + 6];
            const TF r32 = rot[rdim * n + 7];
            const TF r33 = rot[rdim * n + 8];

            TF u_ip = TF(0.);
            TF v_ip = TF(0.);
            TF w_ip = TF(0.);

            // 1. interpolate surroundings neighbours to interpolation point
            for (int i = 0; i < n_idw_loc; ++i)
            {
                const int ii = i + n * n_idw_loc;                                  
                const int ijku = ipui[ii] + ipuj[ii] * icells + ipuk[ii] * ijcells;
                const int ijkv = ipvi[ii] + ipvj[ii] * icells + ipvk[ii] * ijcells;
                const int ijkw = ipwi[ii] + ipwj[ii] * icells + ipwk[ii] * ijcells;

                u_ip += c_idw_u[ii] * fld_u[ijku];
                v_ip += c_idw_v[ii] * fld_v[ijkv];
                w_ip += c_idw_w[ii] * fld_w[ijkw];
            }

            // 2. rotate velocities to locally align with surface tangent (under the assumption that flow at second layer still aligns)
            u_ip_la = r11 * u_ip + r12 * v_ip + r13 * w_ip;
            v_ip_la = r21 * u_ip + r22 * v_ip + r23 * w_ip;

            umag_ip_la = std::pow(fm::pow2(u_ip_la) + fm::pow2(v_ip_la), TF(0.5));

            if (db[n] > z0b[n])
            {
                // 3. calculate local shear velocity (ustar) and thereby eddy viscosity at forcing point
                // for now, (1) neglect flow rotation over height, (2) neglect stability effects (requires "fine enough" grid),
                // (3) assume both points are in logarithmic layer, and (4) assume zero-valued Dirichlet conditions for momentum (i.e. no-slip condition)
                // future options: investigate use of Van Driest (1956) correction and/or DNS mode.
                fld_evisc[ijkf] = fm::pow2(Constants::kappa<TF>) * umag_ip_la * db[n] / std::log(di[n] / z0b[n]);
            }
            else
            {
                // SvdL, 20240901: Currently defaulting to WRONG(?) evisc-value due to zero boundary distance..
                fld_evisc[ijkf] = fm::pow2(Constants::kappa<TF>) * umag_ip_la * z0b[n] / std::log(di[n] / z0b[n]);
            }
        }
    }


    template<typename TF> __global__
    void set_ib_points_g(
        TF* const __restrict__ tend_var,
        const TF* const __restrict__ var,
        const int* const __restrict__ ijk_ib, 
        const TF val,
        const int n_ib,
        const double dt)
    {
        const int n    = blockIdx.x*blockDim.x + threadIdx.x;
        const TF  dtf  = static_cast<TF>(dt); // SvdL, 20240901: seems like ugly solution, see how to improve later.

        if (n < n_ib)
        {
            tend_var[ijk_ib[n]] = ( val - var[ijk_ib[n]] ) / dtf;       
        }

    }

    template<typename TF> __global__
    void set_ib_points_evisc_g(
        TF* const __restrict__ var,
        const int* const __restrict__ ijk_ib, 
        const TF val,
        const int n_ib)
    {
        const int n = blockIdx.x*blockDim.x + threadIdx.x;

        if (n < n_ib)
            var[ijk_ib[n]] = val;
        
    }
}

// SvdL, 20240908: started CUDA implementation
#ifdef USECUDA

template <typename TF>
void Immersed_boundary<TF>::exec_viscosity()
{
    if (sw_ib == IB_type::Disabled)
    return;

    const int blocki = 256;

    const int n_fp_c = fpoints.at("s").n_fpoints;
    const int n_ib_c = fpoints.at("s").n_ibpoints;

    const int gridfp_c = n_fp_c / blocki + (n_fp_c % blocki > 0);
    const int gridib_c = n_ib_c / blocki + (n_ib_c % blocki > 0);

    dim3 gridGPU_fp_c(gridfp_c);
    dim3 gridGPU_ib_c(gridib_c);

    dim3 blockGPU(blocki);

    set_forcing_points_evisc_g<TF><<gridGPU_fp_c, blockGPU>>>(
            fields.sd.at("evisc")->fld_g,
            fields.mp.at("u")->fld_g,
            fields.mp.at("v")->fld_g,
            fields.mp.at("w")->fld_g,
            // fpoints.at("s").sbot.at("s")_g,                                                                       //<< SvdL, 20240909: for now not needed here,  // value of boundary conditions to be enforced
            fpoints.at("s").i_g, fpoints.at("s").j_g, fpoints.at("s").k_g,                                           // points to be forced
            fpoints.at("s").rot_g,                                                                                   // rotational matrix for local surface alignment
            fpoints.at("s").ip_u_i_g, fpoints.at("s").ip_u_j_g, fpoints.at("s").ip_u_k_g, fpoints.at("s").c_idw_u_g, // locations of the neighbouring u-points + weights
            fpoints.at("s").ip_v_i_g, fpoints.at("s").ip_v_j_g, fpoints.at("s").ip_v_k_g, fpoints.at("s").c_idw_v_g, // locations of the neighbouring v-points + weights
            fpoints.at("s").ip_w_i_g, fpoints.at("s").ip_w_j_g, fpoints.at("s").ip_w_k_g, fpoints.at("s").c_idw_w_g, // locations of the neighbouring w-points + weights
            fpoints.at("s").ip_s_i_g, fpoints.at("s").ip_s_j_g, fpoints.at("s").ip_s_k_g, fpoints.at("s").c_idw_s_g, // locations of the neighbouring s-points + weights
            fpoints.at("s").dist_b_g,                                                                                // distance nearest immersed boundary point to forcing point
            fpoints.at("s").dist_i_g,                                                                                // distance interpolation point to forcing point
            fpoints.at("s").z0b_g,                                                                                   // local roughness lengths of forcing points (all scalars will have same for now..)
            Boundary_type::Dirichlet_type,                                                                           // should contain Boundary_Type:: for all scalars (make variation between scalars possible?), also unused for evisc
            fields.visc, n_fp_c, this->n_idw_points,
            gd.icells, gd.ijcells);
    cuda_check_error();

    // eddy viscosity should just be zero inside objects
    set_ib_points_evisc_g<TF><<gridGPU_fib_c, blockGPU>>>(
            fields.sd.at("evisc")->fld_g,
            ibpoints.at("s").ijk._g, TF(0.), 
            n_ib_c);
    cuda_check_error();

    // Enforce cyclic boundary conditions for updated evisc
    boundary_cyclic.exec_g(fields.sd.at("evisc")->fld_g);

}

template <typename TF>
void Immersed_boundary<TF>::exec(const double dt)
{
    if (sw_ib == IB_type::Disabled)
    return;

    const int blocki = 256;

    const int n_fp_u = fpoints.at("u").n_fpoints;
    const int n_ib_u = fpoints.at("u").n_ibpoints;
    const int n_fp_v = fpoints.at("v").n_fpoints;
    const int n_ib_v = fpoints.at("v").n_ibpoints;
    const int n_fp_w = fpoints.at("w").n_fpoints;
    const int n_ib_w = fpoints.at("w").n_ibpoints;
    const int n_fp_c = fpoints.at("s").n_fpoints;
    const int n_ib_c = fpoints.at("s").n_ibpoints;

    const int gridfp_u = n_fp_u / blocki + (n_fp_u % blocki > 0);
    const int gridib_u = n_ib_u / blocki + (n_ib_u % blocki > 0);
    const int gridfp_v = n_fp_v / blocki + (n_fp_v % blocki > 0);
    const int gridib_v = n_ib_v / blocki + (n_ib_v % blocki > 0);
    const int gridfp_w = n_fp_w / blocki + (n_fp_w % blocki > 0);
    const int gridib_w = n_ib_w / blocki + (n_ib_w % blocki > 0);
    const int gridfp_c = n_fp_c / blocki + (n_fp_c % blocki > 0);
    const int gridib_c = n_ib_c / blocki + (n_ib_c % blocki > 0);

    dim3 gridGPU_fp_u(gridfp_u);
    dim3 gridGPU_ib_u(gridib_u);
    dim3 gridGPU_fp_v(gridfp_v);
    dim3 gridGPU_ib_v(gridib_v);
    dim3 gridGPU_fp_w(gridfp_w);       
    dim3 gridGPU_ib_w(gridib_w);
    dim3 gridGPU_fp_c(gridfp_c);
    dim3 gridGPU_ib_c(gridib_c);

    dim3 blockGPU(blocki);

    set_forcing_points_u_g<TF><<<gridGPU_fp_u, blockGPU>>>(
            fields.mt.at("u")->fld_g,
            fields.mt.at("v")->fld_g,
            fields.mt.at("w")->fld_g,
            fields.mp.at("u")->fld_g,
            fields.mp.at("v")->fld_g,
            fields.mp.at("w")->fld_g,
            fpoints.at("u").mbot_g,                                                                                            // value of boundary conditions to be enforced
            fpoints.at("u").i_g, fpoints.at("u").j_g, fpoints.at("u").k_g,                                           // points to be forced
            fpoints.at("u").rot_g,                                                                                             // rotational matrix for local surface alignment
            fpoints.at("u").ip_u_i_g, fpoints.at("u").ip_u_j_g, fpoints.at("u").ip_u_k_g, fpoints.at("u").c_idw_u_g, // locations of the neighbouring u-points + weights
            fpoints.at("u").ip_v_i_g, fpoints.at("u").ip_v_j_g, fpoints.at("u").ip_v_k_g, fpoints.at("u").c_idw_v_g, // locations of the neighbouring v-points + weights
            fpoints.at("u").ip_w_i_g, fpoints.at("u").ip_w_j_g, fpoints.at("u").ip_w_k_g, fpoints.at("u").c_idw_w_g, // locations of the neighbouring w-points + weights
            fpoints.at("u").ip_s_i_g, fpoints.at("u").ip_s_j_g, fpoints.at("u").ip_s_k_g, fpoints.at("u").c_idw_s_g, // locations of the neighbouring s-points + weights
            fpoints.at("u").dist_b_g,                                                                                              // distance nearest immersed boundary point to forcing point
            fpoints.at("u").dist_i_g,                                                                                              // distance interpolation point to forcing point
            fpoints.at("u").z0b_g,                                                                                             // local roughness lengths of forcing points
            Boundary_type::Dirichlet_type,                                                                                          // only allow no-slip conditions for momentum (for now..)
            fields.visc, n_fp_u, this->n_idw_points,
            gd.icells, gd.ijcells,
            dt);
    cuda_check_error();

    set_ib_points_g<TF><<<gridGPU_ib_u, blockGPU>>>(
            fields.mt.at("u")->fld_g,
            fields.mp.at("u")->fld_g,
            ibpoints.at("u").ijk_g, TF(0.), 
            n_ib_u,
            dt);
    cuda_check_error();

    set_forcing_points_v_g<TF><<<gridGPU_fp_v, blockGPU>>>(
            fields.mt.at("u")->fld_g,
            fields.mt.at("v")->fld_g,
            fields.mt.at("w")->fld_g,
            fields.mp.at("u")->fld_g,
            fields.mp.at("v")->fld_g,
            fields.mp.at("w")->fld_g,
            fpoints.at("v").mbot_g,                                                                    // value of boundary conditions to be enforced
            fpoints.at("v").i_g, fpoints.at("v").j_g, fpoints.at("v").k_g,                   // points to be forced
            fpoints.at("v").rot_g,                                                                     // rotational matrix for local surface alignment
            fpoints.at("v").ip_u_i_g, fpoints.at("u").ip_u_j_g, fpoints.at("v").ip_u_k_g, fpoints.at("v").c_idw_u_g, // locations of the neighbouring u-points + weights
            fpoints.at("v").ip_v_i_g, fpoints.at("v").ip_v_j_g, fpoints.at("v").ip_v_k_g, fpoints.at("v").c_idw_v_g, // locations of the neighbouring v-points + weights
            fpoints.at("v").ip_w_i_g, fpoints.at("v").ip_w_j_g, fpoints.at("v").ip_w_k_g, fpoints.at("v").c_idw_w_g, // locations of the neighbouring w-points + weights
            fpoints.at("v").ip_s_i_g, fpoints.at("v").ip_s_j_g, fpoints.at("v").ip_s_k_g, fpoints.at("v").c_idw_s_g, // locations of the neighbouring s-points + weights
            fpoints.at("v").dist_b_g,                                                                      // distance nearest immersed boundary point to forcing point
            fpoints.at("v").dist_i_g,                                                                      // distance interpolation point to forcing point
            fpoints.at("v").z0b_g,                                                                     // local roughness lengths of forcing points
            Boundary_type::Dirichlet_type,                                                                  // only allow no-slip conditions for momentum (for now..)
            fields.visc, n_fp_v, this->n_idw_points,
            gd.icells, gd.ijcells, 
            dt);
    cuda_check_error();

    set_ib_points_g<TF><<<gridGPU_ib_v, blockGPU>>>(
            fields.mt.at("v")->fld_g,
            fields.mp.at("v")->fld_g,
            ibpoints.at("v").ijk_g, TF(0.), 
            n_ib_v,
            dt);
    cuda_check_error();

    set_forcing_points_w_g<TF><<<gridGPU_fp_w, blockGPU>>>(
            fields.mt.at("u")->fld_g,
            fields.mt.at("v")->fld_g,
            fields.mt.at("w")->fld_g,
            fields.mp.at("u")->fld_g,
            fields.mp.at("v")->fld_g,
            fields.mp.at("w")->fld_g,
            fpoints.at("w").mbot_g,                                                                    // value of boundary conditions to be enforced
            fpoints.at("w").i_g, fpoints.at("w").j_g, fpoints.at("w").k_g,                   // points to be forced
            fpoints.at("w").rot_g,                                                                     // rotational matrix for local surface alignment
            fpoints.at("w").ip_u_i_g, fpoints.at("w").ip_u_j_g, fpoints.at("w").ip_u_k_g, fpoints.at("u").c_idw_u_g, // locations of the neighbouring u-points + weights
            fpoints.at("w").ip_v_i_g, fpoints.at("w").ip_v_j_g, fpoints.at("w").ip_v_k_g, fpoints.at("u").c_idw_v_g, // locations of the neighbouring v-points + weights
            fpoints.at("w").ip_w_i_g, fpoints.at("w").ip_w_j_g, fpoints.at("w").ip_w_k_g, fpoints.at("u").c_idw_w_g, // locations of the neighbouring w-points + weights
            fpoints.at("w").ip_s_i_g, fpoints.at("w").ip_s_j_g, fpoints.at("w").ip_s_k_g, fpoints.at("u").c_idw_s_g, // locations of the neighbouring s-points + weights
            fpoints.at("w").dist_b_g,                                                                      // distance nearest immersed boundary point to forcing point
            fpoints.at("w").dist_i_g,                                                                      // distance interpolation point to forcing point
            fpoints.at("w").z0b_g,                                                                     // local roughness lengths of forcing points
            Boundary_type::Dirichlet_type,                                                                  // only allow no-slip conditions for momentum (for now..)
            fields.visc, n_fp_w, this->n_idw_points,
            gd.icells, gd.ijcells,
            dt);
    cuda_check_error();    

    set_ib_points_g<TF><<<gridGPU_ib_w, blockGPU>>>(
            fields.mt.at("w")->fld_g,
            fields.mp.at("w")->fld_g,
            ibpoints.at("w").ijk_g, TF(0.), 
            n_ib_w,
            dt);
    cuda_check_error();

    for (auto &it : fields.sp)
    {
        set_forcing_points_c_g<TF><<<gridGPU_fp_c, blockGPU>>>(
                fields.st.at(it.first)->fld_g,
                it.second->fld_g,
                fields.mp.at("u")->fld_g,
                fields.mp.at("v")->fld_g,
                fields.mp.at("w")->fld_g,
                fpoints.at("s").sbot.at(it.first)_g,                                                       // value of boundary conditions to be enforced
                fpoints.at("s").i_g, fpoints.at("s").j_g, fpoints.at("s").k_g,                   // points to be forced
                fpoints.at("s").rot_g,                                                                     // rotational matrix for local surface alignment (although not used yet for scalars: DO calculate, needed for evisc)
                fpoints.at("s").ip_u_i_g, fpoints.at("s").ip_u_j_g, fpoints.at("s").ip_u_k_g, fpoints.at("s").c_idw_u_g, // locations of the neighbouring u-points + weights
                fpoints.at("s").ip_v_i_g, fpoints.at("s").ip_v_j_g, fpoints.at("s").ip_v_k_g, fpoints.at("s").c_idw_v_g, // locations of the neighbouring v-points + weights
                fpoints.at("s").ip_w_i_g, fpoints.at("s").ip_w_j_g, fpoints.at("s").ip_w_k_g, fpoints.at("s").c_idw_w_g, // locations of the neighbouring w-points + weights
                fpoints.at("s").ip_s_i_g, fpoints.at("s").ip_s_j_g, fpoints.at("s").ip_s_k_g, fpoints.at("s").c_idw_s_g, // locations of the neighbouring s-points + weights
                fpoints.at("s").dist_b_g,                                                                      // distance nearest immersed boundary point to forcing point
                fpoints.at("s").dist_i_g,                                                                      // distance interpolation point to forcing point
                fpoints.at("s").z0b_g,                                                                     // local roughness lengths of forcing points (all scalars will have same for now..)
                sbc.at(it.first),                                                                                         // should contain Boundary_Type:: for all scalars (make variation between scalars possible?)
                fields.visc, n_fp_c, this->n_idw_points,
                gd.icells, gd.ijcells,
                dt);
        cuda_check_error();

        set_ib_points_g<TF><<<gridGPU_ib_c, blockGPU>>>(
                fields.st.at(it.first)->fld_g,
                it.second->fld_g,
                ibpoints.at("s").ijk_g, TF(0.), 
                n_ib_c,
                dt); 
        cuda_check_error(); 
    }  
}
#endif

template <typename TF>
void Immersed_boundary<TF>::prepare_device()
{
    if (sw_ib == IB_type::Disabled)
        return;

    // Allocate and forward FP information (forcing points)
    for (auto& fp : fpoints)
    {
        const int n_fpoints = fp.second.n_fpoints;
        const int rdim = 9;

        const int imemsize_1d = n_fpoints*sizeof(int);
        const int fmemsize_1d = n_fpoints*sizeof(TF);

        const int fmemsize_ro = n_fpoints*rdim*sizeof(TF);

        const int imemsize_2d = n_fpoints*n_idw_points*sizeof(int);
        const int fmemsize_2d = n_fpoints*n_idw_points*sizeof(TF);

         // Allocate
        cuda_safe_call(cudaMalloc(&fp.second.i_g  , imemsize_1d));
        cuda_safe_call(cudaMalloc(&fp.second.j_g  , imemsize_1d));
        cuda_safe_call(cudaMalloc(&fp.second.k_g  , imemsize_1d));
        // cuda_safe_call(cudaMalloc(&fp.second.ijk_g, imemsize_1d));

        cuda_safe_call(cudaMalloc(&fp.second.rot_g   , fmemsize_ro));

        cuda_safe_call(cudaMalloc(&fp.second.dist_b_g, fmemsize_1d));
        cuda_safe_call(cudaMalloc(&fp.second.dist_i_g, fmemsize_1d));
        cuda_safe_call(cudaMalloc(&fp.second.z0b_g   , fmemsize_1d));

        cuda_safe_call(cudaMalloc(&fp.second.ip_u_i_g, imemsize_2d));
        cuda_safe_call(cudaMalloc(&fp.second.ip_u_j_g, imemsize_2d));
        cuda_safe_call(cudaMalloc(&fp.second.ip_u_k_g, imemsize_2d));
        cuda_safe_call(cudaMalloc(&fp.second.ip_v_i_g, imemsize_2d));
        cuda_safe_call(cudaMalloc(&fp.second.ip_v_j_g, imemsize_2d));
        cuda_safe_call(cudaMalloc(&fp.second.ip_v_k_g, imemsize_2d));
        cuda_safe_call(cudaMalloc(&fp.second.ip_w_i_g, imemsize_2d));
        cuda_safe_call(cudaMalloc(&fp.second.ip_w_j_g, imemsize_2d));
        cuda_safe_call(cudaMalloc(&fp.second.ip_w_k_g, imemsize_2d));
        cuda_safe_call(cudaMalloc(&fp.second.ip_s_i_g, imemsize_2d));
        cuda_safe_call(cudaMalloc(&fp.second.ip_s_j_g, imemsize_2d));
        cuda_safe_call(cudaMalloc(&fp.second.ip_s_k_g, imemsize_2d));

        cuda_safe_call(cudaMalloc(&fp.second.c_idw_u_g, fmemsize_2d));
        cuda_safe_call(cudaMalloc(&fp.second.c_idw_v_g, fmemsize_2d));
        cuda_safe_call(cudaMalloc(&fp.second.c_idw_w_g, fmemsize_2d));
        cuda_safe_call(cudaMalloc(&fp.second.c_idw_s_g, fmemsize_2d));

        if (fp.first == "u" || fp.first == "v" || fp.first == "w")
            cuda_safe_call(cudaMalloc(&fp.second.mbot_g, fmemsize_1d));
        else
            for (auto& it : fp.second.sbot)
                cuda_safe_call(cudaMalloc(&fp.second.sbot_g[it.first], fmemsize_1d));

        // Forward copy
        cuda_safe_call(cudaMemcpy(fp.second.i_g, fp.second.i.data(), imemsize_1d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(fp.second.j_g, fp.second.j.data(), imemsize_1d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(fp.second.k_g, fp.second.k.data(), imemsize_1d, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(fp.second.dist_b_g, fp.second.dist_b.data(), fmemsize_1d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(fp.second.dist_i_g, fp.second.dist_i.data(), fmemsize_1d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(fp.second.z0b_g   , fp.second.z0b.data()   , fmemsize_1d, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(fp.second.rot_g   , fp.second.rot.data()   , fmemsize_ro, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(fp.second.ip_u_i_g, fp.second.ip_u_i.data(), imemsize_2d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(fp.second.ip_u_j_g, fp.second.ip_u_j.data(), imemsize_2d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(fp.second.ip_u_k_g, fp.second.ip_u_k.data(), imemsize_2d, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(fp.second.ip_v_i_g, fp.second.ip_v_i.data(), imemsize_2d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(fp.second.ip_v_j_g, fp.second.ip_v_j.data(), imemsize_2d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(fp.second.ip_v_k_g, fp.second.ip_v_k.data(), imemsize_2d, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(fp.second.ip_w_i_g, fp.second.ip_w_i.data(), imemsize_2d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(fp.second.ip_w_j_g, fp.second.ip_w_j.data(), imemsize_2d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(fp.second.ip_w_k_g, fp.second.ip_w_k.data(), imemsize_2d, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(fp.second.ip_s_i_g, fp.second.ip_s_i.data(), imemsize_2d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(fp.second.ip_s_j_g, fp.second.ip_s_j.data(), imemsize_2d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(fp.second.ip_s_k_g, fp.second.ip_s_k.data(), imemsize_2d, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(fp.second.c_idw_u_g, fp.second.c_idw_u.data(), fmemsize_2d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(fp.second.c_idw_v_g, fp.second.c_idw_v.data(), fmemsize_2d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(fp.second.c_idw_w_g, fp.second.c_idw_w.data(), fmemsize_2d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(fp.second.c_idw_s_g, fp.second.c_idw_s.data(), fmemsize_2d, cudaMemcpyHostToDevice));

        if (fp.first == "u" || fp.first == "v" || fp.first == "w")
            cuda_safe_call(cudaMemcpy(fp.second.mbot_g, fp.second.mbot.data(), fmemsize_1d, cudaMemcpyHostToDevice));
        else
            for (auto& it : fp.second.sbot)
                cuda_safe_call(cudaMemcpy(fp.second.sbot_g.at(it.first), fp.second.sbot.at(it.first).data(), 
                               fmemsize_1d, cudaMemcpyHostToDevice));
    }

    // Allocate and forward IB information (internal points)
    for (auto& ib : ibpoints)
    {
        const int n_ibpoints = ib.second.n_ibpoints;

        const int imemsize_1d = n_ibpoints*sizeof(int);

         // Allocate
        cuda_safe_call(cudaMalloc(&ib.second.i_g  , imemsize_1d));
        cuda_safe_call(cudaMalloc(&ib.second.j_g  , imemsize_1d));
        cuda_safe_call(cudaMalloc(&ib.second.k_g  , imemsize_1d));
        cuda_safe_call(cudaMalloc(&ib.second.ijk_g, imemsize_1d));

        // Forward copy
        cuda_safe_call(cudaMemcpy(ib.second.i_g  , ib.second.i.data()  , imemsize_1d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(ib.second.j_g  , ib.second.j.data()  , imemsize_1d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(ib.second.k_g  , ib.second.k.data()  , imemsize_1d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(ib.second.ijk_g, ib.second.ijk.data(), imemsize_1d, cudaMemcpyHostToDevice));
    }      
}

template <typename TF>
void Immersed_boundary<TF>::clear_device()
{
    if (sw_ib == IB_type::Disabled)
        return;

    // De-allocate all FP information (forcing points)
    for (auto& fp : fpoints)
    {
        cuda_safe_call(cudaFree(fp.second.i_g));
        cuda_safe_call(cudaFree(fp.second.j_g));
        cuda_safe_call(cudaFree(fp.second.k_g));

        cuda_safe_call(cudaFree(fp.second.rot_g));

        cuda_safe_call(cudaFree(fp.second.dist_b_g));
        cuda_safe_call(cudaFree(fp.second.dist_i_g));
        cuda_safe_call(cudaFree(fp.second.z0b_g   ));

        cuda_safe_call(cudaFree(fp.second.ip_u_i_g));
        cuda_safe_call(cudaFree(fp.second.ip_u_j_g));
        cuda_safe_call(cudaFree(fp.second.ip_u_k_g));
        cuda_safe_call(cudaFree(fp.second.ip_v_i_g));
        cuda_safe_call(cudaFree(fp.second.ip_v_j_g));
        cuda_safe_call(cudaFree(fp.second.ip_v_k_g));
        cuda_safe_call(cudaFree(fp.second.ip_w_i_g));
        cuda_safe_call(cudaFree(fp.second.ip_w_j_g));
        cuda_safe_call(cudaFree(fp.second.ip_w_k_g));
        cuda_safe_call(cudaFree(fp.second.ip_s_i_g));
        cuda_safe_call(cudaFree(fp.second.ip_s_j_g));
        cuda_safe_call(cudaFree(fp.second.ip_s_k_g));

        cuda_safe_call(cudaFree(fp.second.c_idw_u_g));
        cuda_safe_call(cudaFree(fp.second.c_idw_v_g));
        cuda_safe_call(cudaFree(fp.second.c_idw_w_g));
        cuda_safe_call(cudaFree(fp.second.c_idw_s_g));

        if (fp.first == "u" || fp.first == "v" || fp.first == "w")
            cuda_safe_call(cudaFree(fp.second.mbot_g));
        else
            for (auto& it : fp.second.sbot)
                cuda_safe_call(cudaFree(fp.second.sbot_g[it.first]));
    }

    // De-allocate all IB information (internal points)
    for (auto& ib : ibpoints)
    {
        cuda_safe_call(cudaFree(ib.second.i_g  ));
        cuda_safe_call(cudaFree(ib.second.j_g  ));
        cuda_safe_call(cudaFree(ib.second.k_g  ));
        cuda_safe_call(cudaFree(ib.second.ijk_g));
    }     

}

template class Immersed_boundary<double>;
template class Immersed_boundary<float>;
