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

#include <unistd.h> //<< for DEBUG, SVDL
#include <iostream>
#include <cmath>
#include <algorithm>

#include <constants.h>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "input.h"
#include "immersed_boundary.h"
#include "fast_math.h"
#include "stats.h"
#include "netcdf_interface.h"
#include "cross.h"
#include "finite_difference.h"

using namespace Finite_difference::O2;
namespace
{
    namespace fm = Fast_math;
    
    template<typename TF>
    TF absolute_distance(
            const TF x1, const TF y1, const TF z1,
            const TF x2, const TF y2, const TF z2)
    {
        return std::pow(Fast_math::pow2(x2-x1) + Fast_math::pow2(y2-y1) + Fast_math::pow2(z2-z1), TF(0.5));
    }

    // Help function for sorting std::vector with Neighbour points
    template<typename TF>
    bool compare_value(const Neighbour<TF>& a, const Neighbour<TF>& b)
    {
        return a.distance < b.distance;
    }

    // // SvdL, 20240918: commented out for now, original use unclear
    // bool has_ending(const std::string& full_string, const std::string& ending)
    // {
    //     if (full_string.length() >= ending.length())
    //         return (0 == full_string.compare(full_string.length() - ending.length(), ending.length(), ending));
    //     else
    //         return false;
    // };
        
    template <typename TF>
    void setup_interpolation(
        std::vector<TF>& xi_fp, std::vector<TF>& yi_fp, std::vector<TF>& zi_fp,                             // location of interpolation point
        std::vector<int>& ipui, std::vector<int>& ipuj, std::vector<int>& ipuk, std::vector<TF>& c_idw_u,   // indices of interpolation locations + weights
        std::vector<int>& ipvi, std::vector<int>& ipvj, std::vector<int>& ipvk, std::vector<TF>& c_idw_v, 
        std::vector<int>& ipwi, std::vector<int>& ipwj, std::vector<int>& ipwk, std::vector<TF>& c_idw_w,
        std::vector<int>& ipsi, std::vector<int>& ipsj, std::vector<int>& ipsk, std::vector<TF>& c_idw_s,  
        std::vector<int>& i_fp, std::vector<int>& j_fp, std::vector<int>& k_fp,                             // grid locations of points to be forced
        const int n_fpoints, const int n_idw_points, const int n_idw_points_min,
        TF* const restrict sdfu, TF* const restrict sdfv, TF* const restrict sdfw, TF* const restrict sdfs, // SDFs for all grid locations
        const std::vector<TF>& x, const std::vector<TF>& y, const std::vector<TF>& z, 
        const std::vector<TF>& xh, const std::vector<TF>& yh, const std::vector<TF>& zh,
        const TF dx, const TF dy, 
        const std::vector<TF>& dz, const std::vector<TF>& dzh, 
        const int istart, const int jstart, const int kstart,
        const int iend,   const int jend,   const int kend,
        const int icells, const int ijcells)
    {

        TF c_idw_sum = TF(0.);

        for (int nn=0; nn<n_fpoints; ++nn)
        {
            const TF xi = xi_fp[nn];
            const TF yi = yi_fp[nn];
            const TF zi = zi_fp[nn];

            const int kfp = k_fp[nn];

            // Calculate Delta as "originally" used at forcing point location (should be fine as long as stretching is not too strong)
            const TF Delta  = std::pow(dx*dy*dzh[kfp], TF(1./3.));  // for center positions (mind use of dzh for center position -> see grid definitions)
            const TF Deltah = std::pow(dx*dy* dz[kfp], TF(1./3.));  // for face positions (idem here)

            // Maximum distance allowed for use in inverse distance interpolation (aim to keep interpolation as local as possible)
            const TF dist_max  = TF(2.)*Delta;  //<< value of 2 should be well within range of identified potential neighbours (see below)
            const TF dist_maxh = TF(2.)*Deltah; //<< value of 2 should be well within range of identified potential neighbours (see below)

            // Vectors including all neighbours around interpolation point
            std::vector<Neighbour<TF>> u_neighbours; // neighbouring u-momentum grid points of current (xi,yi,zi) point
            std::vector<Neighbour<TF>> v_neighbours; // idem for v
            std::vector<Neighbour<TF>> w_neighbours; // idem for w
            std::vector<Neighbour<TF>> s_neighbours; // idem for scalars

            // Calculate "semi-nearest" grid indices (i,j) of interpolation point (xi,yi,zi) 
            const int in  = static_cast<int>(std::round( (xi-x[istart]) / dx) ) + istart; // guess of center index i
            const int jn  = static_cast<int>(std::round( (yi-y[istart]) / dy) ) + jstart; // guess of center index j
           
            const int inh = static_cast<int>(std::round( (xi-xh[istart]) / dx) ) + istart; // guess of face index i (used for u points)
            const int jnh = static_cast<int>(std::round( (yi-yh[istart]) / dy) ) + jstart; // guess of face index j (used for v points)

            // Do slightly different procedure for k (in case of stretched grid)
            int kn  = kstart;   // guess of center index k
            int knh = kstart;   // guess of face index k (used for w points)

            for (int k=kstart; k<kend-1; ++k)
            {
                if ( z[k] >= zi ) 
                {
                    if ((z[k]-zi) < (z[k+1]-zi))
                    {
                        kn = k; 
                    }
                    else
                    {
                        kn = k+1; 
                    }
                    break;
                }
            }

            for (int k=kstart; k<kend; ++k)
            {
                if ( zh[k] >= zi ) 
                {
                    if ((zh[k]-zi) < (zh[k+1]-zi))
                    {
                        knh = k; 
                    }
                    else
                    {
                        knh = k+1; 
                    }
                    break;
                }
            }

            // Limit vertical stencil near surface (interpolation locations may not be below surface)
            const int dk0  = std::max(-2, kstart-kn);
            const int dk0h = std::max(-2, kstart-knh);

            // Search for potential neighbours (at u, v, w and center positions) for all forcing points
            // --> availability of SDF enables to explicitly test if potential neighbour is (another) forcing cell to be skipped
            // --> also, enables to preselect points that are not too far from forcing point (e.g., within 4*Delta)

            // 1. Do for uloc

            c_idw_sum = 0.; // reset to zero

            // Find neighbouring grid points outside IB (use rather wide search region)
            for (int dk=dk0; dk<3; ++dk)
                for (int dj=-2; dj<3; ++dj)
                    for (int di=-2; di<3; ++di)
                    {
                        const int ijk_test = (inh+di) + (jn+dj)*icells + (kn+dk)*ijcells; // combined grid index to test
                        
                        if (sdfu[ijk_test] < Delta || sdfu[ijk_test] > TF(4.)*Delta)
                            continue;

                        const TF distance = absolute_distance(xi, yi, zi, xh[inh+di], y[jn+dj], z[kn+dk]);
                        Neighbour<TF> tmp_neighbour = {inh+di, jn+dj, kn+dk, distance};
                        u_neighbours.push_back(tmp_neighbour);
                        
                    }

            // Sort them on distance
            std::sort(u_neighbours.begin(), u_neighbours.end(), compare_value<TF>);

            // If smallest distance is zero (to within some precision), this point gets weight 1 and the rest is set to zero weight.
            if (u_neighbours[0].distance < TF(1e-7))
            {
                ipui[nn*n_idw_points]    = u_neighbours[0].i;
                ipuj[nn*n_idw_points]    = u_neighbours[0].j;
                ipuk[nn*n_idw_points]    = u_neighbours[0].k;
                c_idw_u[nn*n_idw_points] = TF(1.);  

                for (int ii=1; ii<n_idw_points; ++ii)
                {
                    const int in = ii + nn*n_idw_points; // index gives position in vectors

                    ipui[in]    = inh;      // just a "FillValue", make sure it carries weight ZERO
                    ipuj[in]    = jn;
                    ipuk[in]    = kn;
                    c_idw_u[in] = TF(0.);   // set weights here to ZERO
                }
            }
            else
            {
                if (u_neighbours.size() < n_idw_points_min)
                {
                    std::cout << "ERROR: only found " << u_neighbours.size() << " u interpolation points around (less than minimum) ";
                    std::cout << "xi=" << xi << ", yi=" << yi << ", zi=" << zi << std::endl; 
                    throw 1;
                }
                else if (u_neighbours.size() < n_idw_points)
                {
                    std::cout << "NOTE: only found " << u_neighbours.size() << " u interpolation points around (less than maximum allowed) ";
                    std::cout << "xi=" << xi << ", yi=" << yi << ", zi=" << zi << std::endl; 
                }

                for (int ii=0; ii<n_idw_points; ++ii)
                {
                    const int in = ii + nn*n_idw_points; // index represent number of Forcing Points for this variable (gives position in vector)
                    if (ii < u_neighbours.size())
                    {   
                        ipui[in]    = u_neighbours[ii].i;
                        ipuj[in]    = u_neighbours[ii].j;
                        ipuk[in]    = u_neighbours[ii].k;

                        // SvdL, 20240902: using modified Shepard's Method for weights, with p=0.5 (appropriate value?)
                        c_idw_u[in] = std::pow( std::max(TF(0.), (dist_max - u_neighbours[ii].distance)) / (dist_max * u_neighbours[ii].distance) , TF(0.5) );
                        c_idw_sum  += c_idw_u[in];
                    }
                    else
                    {
                        ipui[in]    = inh;      // just a "FillValue", make sure it carries weight ZERO
                        ipuj[in]    = jn;
                        ipuk[in]    = kn;
                        c_idw_u[in] = TF(0.);   // set weights here to ZERO
                    }
                }
                
                // normalize all weights
                for (int ii=0; ii<n_idw_points; ++ii)
                {
                    const int in = ii + nn*n_idw_points; 
                    c_idw_u[in] /= c_idw_sum; 
                }
            }

            // 2. Do for vloc

            c_idw_sum = 0.; // reset to zero

            // Find neighbouring grid points outside IB (use rather wide search region)
            for (int dk=dk0; dk<3; ++dk)
                for (int dj=-2; dj<3; ++dj)
                    for (int di=-2; di<3; ++di)
                    {
                        const int ijk_test = (in+di) + (jnh+dj)*icells + (kn+dk)*ijcells; // combined grid index to test
                        
                        if (sdfv[ijk_test] < Delta || sdfv[ijk_test] > TF(4.)*Delta)
                            continue;

                        const TF distance = absolute_distance(xi, yi, zi, x[in+di], yh[jnh+dj], z[kn+dk]);
                        Neighbour<TF> tmp_neighbour = {in+di, jnh+dj, kn+dk, distance};
                        v_neighbours.push_back(tmp_neighbour);
                        
                    }

            // Sort them on distance
            std::sort(v_neighbours.begin(), v_neighbours.end(), compare_value<TF>);

            // If smallest distance is zero (to within some precision), this point gets weight 1 and the rest is set to zero weight.
            if (v_neighbours[0].distance < TF(1e-7))
            {
                ipvi[nn*n_idw_points]    = v_neighbours[0].i;
                ipvj[nn*n_idw_points]    = v_neighbours[0].j;
                ipvk[nn*n_idw_points]    = v_neighbours[0].k;
                c_idw_v[nn*n_idw_points] = TF(1.);  

                for (int ii=1; ii<n_idw_points; ++ii)
                {
                    const int in = ii + nn*n_idw_points; // index gives position in vectors

                    ipvi[in]    = in;       // just a "FillValue", make sure it carries weight ZERO
                    ipvj[in]    = jnh;
                    ipvk[in]    = kn;
                    c_idw_v[in] = TF(0.);   // set weights here to ZERO
                }
            }
            else
            {
                if (v_neighbours.size() < n_idw_points_min)
                {
                    std::cout << "ERROR: only found " << v_neighbours.size() << " v interpolation points around (less than minimum) ";
                    std::cout << "xi=" << xi << ", yi=" << yi << ", zi=" << zi << std::endl; 
                    throw 1;
                }
                else if (v_neighbours.size() < n_idw_points)
                {
                    std::cout << "NOTE: only found " << v_neighbours.size() << " v interpolation points around (less than maximum allowed) ";
                    std::cout << "xi=" << xi << ", yi=" << yi << ", zi=" << zi << std::endl; 
                }
                
                for (int ii=0; ii<n_idw_points; ++ii)
                {
                    const int in = ii + nn*n_idw_points; // index represent number of Forcing Points for this variable (gives position in vector)
                    if (ii < v_neighbours.size())
                    {   
                        ipvi[in]    = v_neighbours[ii].i;
                        ipvj[in]    = v_neighbours[ii].j;
                        ipvk[in]    = v_neighbours[ii].k;

                        // SvdL, 20240902: using modified Shepard's Method for weights, with p=0.5 (appropriate value?)
                        c_idw_v[in] = std::pow( std::max(TF(0.), (dist_max - v_neighbours[ii].distance)) / (dist_max * v_neighbours[ii].distance) , TF(0.5) );
                        c_idw_sum  += c_idw_v[in];
                    }
                    else
                    {
                        ipvi[in]    = in;       // just a "FillValue", make sure it carries weight ZERO
                        ipvj[in]    = jnh;
                        ipvk[in]    = kn;
                        c_idw_v[in] = TF(0.);   // set weights here to ZERO
                    }
                }
                
                // normalize all weights
                for (int ii=0; ii<n_idw_points; ++ii)
                {
                    const int in = ii + nn*n_idw_points; 
                    c_idw_v[in] /= c_idw_sum; 
                }
            }

            // 3. Do for wloc

            c_idw_sum = 0.; // reset to zero

            // Find neighbouring grid points outside IB (use rather wide search region)
            for (int dk=dk0h; dk<3; ++dk)
                for (int dj=-2; dj<3; ++dj)
                    for (int di=-2; di<3; ++di)
                    {
                        const int ijk_test = (in+di) + (jn+dj)*icells + (knh+dk)*ijcells; // combined grid index to test
                        
                        if (sdfw[ijk_test] < Deltah || sdfw[ijk_test] > TF(4.)*Deltah)
                            continue;

                        const TF distance = absolute_distance(xi, yi, zi, x[in+di], y[jn+dj], zh[knh+dk]);
                        Neighbour<TF> tmp_neighbour = {in+di, jn+dj, knh+dk, distance};
                        w_neighbours.push_back(tmp_neighbour);
                        
                    }

            // Sort them on distance
            std::sort(w_neighbours.begin(), w_neighbours.end(), compare_value<TF>);

            // If smallest distance is zero (to within some precision), this point gets weight 1 and the rest is set to zero weight.
            if (w_neighbours[0].distance < TF(1e-7))
            {

                ipwi[nn*n_idw_points]    = w_neighbours[0].i;
                ipwj[nn*n_idw_points]    = w_neighbours[0].j;
                ipwk[nn*n_idw_points]    = w_neighbours[0].k;
                c_idw_w[nn*n_idw_points] = TF(1.);  

                for (int ii=1; ii<n_idw_points; ++ii)
                {
                    const int in = ii + nn*n_idw_points; // index gives position in vectors

                    ipwi[in]    = in;       // just a "FillValue", make sure it carries weight ZERO
                    ipwj[in]    = jn;
                    ipwk[in]    = knh;
                    c_idw_w[in] = TF(0.);   // set weights here to ZERO
                }
            }
            else
            {
                if (w_neighbours.size() < n_idw_points_min)
                {
                    std::cout << "ERROR: only found " << w_neighbours.size() << " w interpolation points around (less than minimum) ";
                    std::cout << "xi=" << xi << ", yi=" << yi << ", zi=" << zi << std::endl; 
                    throw 1;
                }
                else if (w_neighbours.size() < n_idw_points)
                {
                    std::cout << "NOTE: only found " << w_neighbours.size() << " w interpolation points around (less than maximum allowed) ";
                    std::cout << "xi=" << xi << ", yi=" << yi << ", zi=" << zi << std::endl; 
                }
                
                for (int ii=0; ii<n_idw_points; ++ii)
                {
                    const int in = ii + nn*n_idw_points; // index represent number of Forcing Points for this variable (gives position in vector)
                    if (ii < w_neighbours.size())
                    {   
                        ipwi[in]    = w_neighbours[ii].i;
                        ipwj[in]    = w_neighbours[ii].j;
                        ipwk[in]    = w_neighbours[ii].k;

                        // SvdL, 20240902: using modified Shepard's Method for weights, with p=0.5 (appropriate value?)
                        c_idw_w[in] = std::pow( std::max(TF(0.), (dist_maxh - w_neighbours[ii].distance)) / (dist_max * w_neighbours[ii].distance) , TF(0.5) );
                        c_idw_sum  += c_idw_w[in];
                    }
                    else
                    {
                        ipui[in]    = in;       // just a "FillValue", make sure it carries weight ZERO
                        ipuj[in]    = jn;
                        ipuk[in]    = knh;
                        c_idw_w[in] = TF(0.);   // set weights here to ZERO
                    }
                }
                
                // normalize all weights
                for (int ii=0; ii<n_idw_points; ++ii)
                {
                    const int in = ii + nn*n_idw_points; 
                    c_idw_w[in] /= c_idw_sum; 
                }
            }

            // 4. Do for sloc

            c_idw_sum = 0.; // reset to zero

            // Find neighbouring grid points outside IB (use rather wide search region)
            for (int dk=dk0; dk<3; ++dk)
                for (int dj=-2; dj<3; ++dj)
                    for (int di=-2; di<3; ++di)
                    {
                        const int ijk_test = (in+di) + (jn+dj)*icells + (kn+dk)*ijcells; // combined grid index to test
                        
                        if (sdfs[ijk_test] < Delta || sdfs[ijk_test] > TF(4.)*Delta)
                            continue;

                        const TF distance = absolute_distance(xi, yi, zi, x[in+di], y[jn+dj], z[kn+dk]);
                        Neighbour<TF> tmp_neighbour = {in+di, jn+dj, kn+dk, distance};
                        s_neighbours.push_back(tmp_neighbour);
                        
                    }

            // Sort them on distance
            std::sort(s_neighbours.begin(), s_neighbours.end(), compare_value<TF>);

            // If smallest distance is zero (to within some precision), this point gets weight 1 and the rest is set to zero weight.
            if (s_neighbours[0].distance < TF(1e-7))
            {

                ipsi[nn*n_idw_points]    = s_neighbours[0].i;
                ipsj[nn*n_idw_points]    = s_neighbours[0].j;
                ipsk[nn*n_idw_points]    = s_neighbours[0].k;
                c_idw_s[nn*n_idw_points] = TF(1.);  

                for (int ii=1; ii<n_idw_points; ++ii)
                {
                    const int in = ii + nn*n_idw_points; // index gives position in vectors

                    ipsi[in]    = in;       // just a "FillValue", make sure it carries weight ZERO
                    ipsj[in]    = jn;
                    ipsk[in]    = kn;
                    c_idw_s[in] = TF(0.);   // set weights here to ZERO
                }
            }
            else
            {
                if (s_neighbours.size() < n_idw_points_min)
                {
                    std::cout << "ERROR: only found " << s_neighbours.size() << " s interpolation points around (less than minimum) ";
                    std::cout << "xi=" << xi << ", yi=" << yi << ", zi=" << zi << std::endl; 
                    throw 1;
                }
                else if (s_neighbours.size() < n_idw_points)
                {
                    std::cout << "NOTE: only found " << s_neighbours.size() << " s interpolation points around (less than maximum allowed) ";
                    std::cout << "xi=" << xi << ", yi=" << yi << ", zi=" << zi << std::endl; 
                }
                
                for (int ii=0; ii<n_idw_points; ++ii)
                {
                    const int in = ii + nn*n_idw_points; // index represent number of Forcing Points for this variable (gives position in vector)
                    if (ii < s_neighbours.size())
                    {   
                        ipsi[in]    = s_neighbours[ii].i;
                        ipsj[in]    = s_neighbours[ii].j;
                        ipsk[in]    = s_neighbours[ii].k;

                        // SvdL, 20240902: using modified Shepard's Method for weights, with p=0.5 (appropriate value?)
                        c_idw_s[in] = std::pow( std::max(TF(0.), (dist_max - s_neighbours[ii].distance)) / (dist_max * s_neighbours[ii].distance) , TF(0.5) );
                        c_idw_sum  += c_idw_s[in];
                    }
                    else
                    {
                        ipsi[in]    = in;       // just a "FillValue", make sure it carries weight ZERO
                        ipsj[in]    = jn;
                        ipsk[in]    = kn;
                        c_idw_s[in] = TF(0.);   // set weights here to ZERO
                    }
                }
                
                // normalize all weights
                for (int ii=0; ii<n_idw_points; ++ii)
                {
                    const int in = ii + nn*n_idw_points; 
                    c_idw_s[in] /= c_idw_sum; 
                }
            }

        }

    }

    template<typename TF>
    void process_sdf(
        Forcing_points<TF>& fpoints, IB_points<TF>& ibpoints, 
        TF* const restrict sdf,
        const std::vector<TF>& x, const std::vector<TF>& y, const std::vector<TF>& z,
        const TF dx, const TF dy, 
        const std::vector<TF>& dz, const std::vector<TF>& dzi,
        const int istart, const int jstart, const int kstart,
        const int iend,   const int jend,   const int kend,
        const int icells, const int jcells, const int ijcells,
        const int n_idw_points)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;
        const TF dxi = TF(1.)/dx;
        const TF dyi = TF(1.)/dy;

        TF dsdx;
        TF dsdy;
        TF dsdz;
        TF norm;
        TF faci;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        const TF Delta = std::pow(dx*dy*dz[k], TF(1./3.)); // Only points in shell distance=Delta are to be forced

                        if ( sdf[ijk] < TF(0.))
                        {
                            ibpoints.i.push_back(i);
                            ibpoints.j.push_back(j);
                            ibpoints.k.push_back(k);
                            ibpoints.ijk.push_back(ijk);
                        } 
                        else if ( TF(0.) <= sdf[ijk] < Delta) // SvdL, 20240828: for now, points on boundary must count as forcing point... points at Delta are "free". 
                        // this may create a problem: 2 forcing points per wall location may be found (if one is on the wall, another one above may be identified)... THINK about this later.
                        {
                            fpoints.i.push_back(i);
                            fpoints.j.push_back(j);
                            fpoints.k.push_back(k);
                            fpoints.ijk.push_back(ijk);
                            fpoints.dist_b.push_back(sdf[ijk]);

                            // Obtain local gradient to find interpolation point and surface point (also equals local normal vector; should point away from surface)
                            dsdx = ( interp2(sdf[ijk], sdf[ijk+ii      ]) - interp2(sdf[ijk-ii      ], sdf[ijk]) ) * dxi;
                            dsdy = ( interp2(sdf[ijk], sdf[ijk   +jj   ]) - interp2(sdf[ijk   -jj   ], sdf[ijk]) ) * dyi; 
                            dsdz = ( interp2(sdf[ijk], sdf[ijk      +kk]) - interp2(sdf[ijk      -kk], sdf[ijk]) ) * dzi[k]; 

                            norm = std::sqrt( fm::pow2(dsdx) + fm::pow2(dsdy) + fm::pow2(dsdz) );

                            dsdx /= norm;
                            dsdy /= norm;
                            dsdz /= norm;

                            // Calculate and store interpolation point (xi,..) and surface point (xb,..)
                            faci = ( sdf[ijk] < TF(0.5)*Delta ) ? (1.25*Delta/sdf[ijk]) : (1.75*Delta/sdf[ijk]) ; //<< This should extend the interpolation point to well within fluid zone (thus surrounded by free fluid points)

                            fpoints.xi.push_back(x[i] + faci*sdf[ijk]*dsdx);
                            fpoints.yi.push_back(y[j] + faci*sdf[ijk]*dsdy);
                            fpoints.zi.push_back(z[k] + faci*sdf[ijk]*dsdz);
                            fpoints.xb.push_back(x[i] - sdf[ijk]*dsdx);
                            fpoints.yb.push_back(y[j] - sdf[ijk]*dsdy);
                            fpoints.zb.push_back(z[k] - sdf[ijk]*dsdz);
                            fpoints.dist_i.push_back((faci+TF(1.))*sdf[ijk]);         //<< i.c.w. calculation of faci, this is ugly but clear

                            // Calculate and store elememts of rotation matrix (as normal vector is now readily available)
                            // ... this whole part may seem sloppy: always possibility to change and calculate inside function where points are forced
                            if (dsdx == TF(1.) && dsdy == TF(0.) && dsdz == TF(0.))
                            {
                                fpoints.rot.push_back(TF(0.));                              
                                fpoints.rot.push_back(TF(0.));
                                fpoints.rot.push_back(TF(-1.));
                                fpoints.rot.push_back(TF(0.));
                                fpoints.rot.push_back(TF(1.));
                                fpoints.rot.push_back(TF(0.));
                                fpoints.rot.push_back(TF(1.));
                                fpoints.rot.push_back(TF(0.));
                                fpoints.rot.push_back(TF(0.));
                            }
                            else if (dsdx == TF(-1.) && dsdy == TF(0.) && dsdz == TF(0.))
                            {
                                fpoints.rot.push_back(TF(0.));                              
                                fpoints.rot.push_back(TF(0.));
                                fpoints.rot.push_back(TF(1.));
                                fpoints.rot.push_back(TF(0.));
                                fpoints.rot.push_back(TF(1.));
                                fpoints.rot.push_back(TF(0.));
                                fpoints.rot.push_back(TF(-1.));
                                fpoints.rot.push_back(TF(0.));
                                fpoints.rot.push_back(TF(0.));
                            }
                            else
                            {
                                const TF scale = std::sqrt( fm::pow2(dsdy) + fm::pow2(dsdz)); 

                                fpoints.rot.push_back((fm::pow2(dsdy) + fm::pow2(dsdz)) / scale);                              
                                fpoints.rot.push_back(-dsdx*dsdy/scale);
                                fpoints.rot.push_back(-dsdx*dsdz/scale);
                                fpoints.rot.push_back(TF(0.));
                                fpoints.rot.push_back(dsdz/scale);
                                fpoints.rot.push_back(-dsdy/scale);
                                fpoints.rot.push_back(dsdx);
                                fpoints.rot.push_back(dsdy);
                                fpoints.rot.push_back(dsdz);
                            }

                        }
                    }

        const int nfpoints = fpoints.i.size();
        fpoints.n_fpoints = nfpoints;
    
        const int nibpoints = ibpoints.i.size();
        ibpoints.n_ibpoints = nibpoints;   

        // Finally, preset containers for SvdL,interpolation points and weights
        // SvdL, 20240902: if we already fill here with "semi"-sensical values, we could potentially reduce unneccesary lines of code in the setup_interpolation function
        fpoints.ip_u_i.resize(nfpoints*n_idw_points);
        fpoints.ip_u_j.resize(nfpoints*n_idw_points);
        fpoints.ip_u_k.resize(nfpoints*n_idw_points);
        fpoints.ip_v_i.resize(nfpoints*n_idw_points);
        fpoints.ip_v_j.resize(nfpoints*n_idw_points);
        fpoints.ip_v_k.resize(nfpoints*n_idw_points);
        fpoints.ip_w_i.resize(nfpoints*n_idw_points);
        fpoints.ip_w_j.resize(nfpoints*n_idw_points);
        fpoints.ip_w_k.resize(nfpoints*n_idw_points);
        fpoints.c_idw_u.resize(nfpoints*n_idw_points);
        fpoints.c_idw_v.resize(nfpoints*n_idw_points);
        fpoints.c_idw_w.resize(nfpoints*n_idw_points);
        fpoints.c_idw_s.resize(nfpoints*n_idw_points);

    }

    template <typename TF>
    void set_forcing_points_u(
        TF* const restrict tend_u,
        const TF* const restrict tend_v,
        const TF* const restrict tend_w,
        const TF* const restrict fld_u,
        const TF* const restrict fld_v,
        const TF* const restrict fld_w,
        const TF* const restrict boundary_value,
        const int* const restrict gi, const int* const restrict gj, const int* const restrict gk,
        const TF* const restrict rot,
        const int* const restrict ipui, const int* const restrict ipuj, const int* const restrict ipuk, const TF* const restrict c_idw_u,
        const int* const restrict ipvi, const int* const restrict ipvj, const int* const restrict ipvk, const TF* const restrict c_idw_v, 
        const int* const restrict ipwi, const int* const restrict ipwj, const int* const restrict ipwk, const TF* const restrict c_idw_w,
        const int* const restrict ipsi, const int* const restrict ipsj, const int* const restrict ipsk, const TF* const restrict c_idw_s, // SvdL, 20240901: not used for now..
        const TF* const restrict db, const TF* const restrict di, const TF* const restrict z0b,
        Boundary_type bc, const TF visc, const int n_fpoints, const int n_idw_loc,
        const int icells, const int ijcells,
        const double dt)
    {

        const int rdim = 9;                                        
        const TF  dtf  = static_cast<TF>(dt); // SvdL, 20240901: seems like ugly solution, see how to improve later.

        TF u_ip_la;
        TF v_ip_la;
        TF w_ip_la;
        TF u_fp_la;
        TF v_fp_la;
        TF w_fp_la;

        // Loop over all points to be forced
        for (int n = 0; n < n_fpoints; ++n)
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
            else // SvdL, 20240901: investigate if can be changed into Van Driest like correction.. would require iterative approach?
            {
                tend_u[ijkf] = ( TF(0.) - fld_u[ijkf] ) / dtf;
            }

        }
    }

    template <typename TF>
    void set_forcing_points_v(
        const TF* const restrict tend_u,
        TF* const restrict tend_v,
        const TF* const restrict tend_w,
        const TF* const restrict fld_u,
        const TF* const restrict fld_v,
        const TF* const restrict fld_w,
        const TF* const restrict boundary_value,
        const int* const restrict gi, const int* const restrict gj, const int* const restrict gk,
        const TF* const restrict rot,
        const int* const restrict ipui, const int* const restrict ipuj, const int* const restrict ipuk, const TF* const restrict c_idw_u,
        const int* const restrict ipvi, const int* const restrict ipvj, const int* const restrict ipvk, const TF* const restrict c_idw_v, 
        const int* const restrict ipwi, const int* const restrict ipwj, const int* const restrict ipwk, const TF* const restrict c_idw_w,
        const int* const restrict ipsi, const int* const restrict ipsj, const int* const restrict ipsk, const TF* const restrict c_idw_s, // SvdL, 20240901: not used for now..
        const TF* const restrict db, const TF* const restrict di, const TF* const restrict z0b,
        Boundary_type bc, const TF visc, const int n_fpoints, const int n_idw_loc,
        const int icells, const int ijcells,
        const double dt)
    {

        const int rdim = 9;                                       
        const TF  dtf  = static_cast<TF>(dt); // SvdL, 20240901: seems like ugly solution, see how to improve later.

        TF u_ip_la;
        TF v_ip_la;
        TF w_ip_la;
        TF u_fp_la;
        TF v_fp_la;
        TF w_fp_la;

        // Loop over all points to be forced
        for (int n = 0; n < n_fpoints; ++n)
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

    template <typename TF>
    void set_forcing_points_w(
        const TF* const restrict tend_u,
        const TF* const restrict tend_v,
        TF* const restrict tend_w,
        const TF* const restrict fld_u,
        const TF* const restrict fld_v,
        const TF* const restrict fld_w,
        const TF* const restrict boundary_value,
        const int* const restrict gi, const int* const restrict gj, const int* const restrict gk,
        const TF* const restrict rot,
        const int* const restrict ipui, const int* const restrict ipuj, const int* const restrict ipuk, const TF* const restrict c_idw_u,
        const int* const restrict ipvi, const int* const restrict ipvj, const int* const restrict ipvk, const TF* const restrict c_idw_v, 
        const int* const restrict ipwi, const int* const restrict ipwj, const int* const restrict ipwk, const TF* const restrict c_idw_w,
        const int* const restrict ipsi, const int* const restrict ipsj, const int* const restrict ipsk, const TF* const restrict c_idw_s, // SvdL, 20240901: not used for now..
        const TF* const restrict db, const TF* const restrict di, const TF* const restrict z0b,
        Boundary_type bc, const TF visc, const int n_fpoints, const int n_idw_loc,
        const int icells, const int ijcells,
        const double dt)
    {

        const int rdim = 9;                                        
        const TF  dtf  = static_cast<TF>(dt); // SvdL, 20240901: seems like ugly solution, see how to improve later.

        TF u_ip_la;
        TF v_ip_la;
        TF w_ip_la;
        TF u_fp_la;
        TF v_fp_la;
        TF w_fp_la;

        // Loop over all points to be forced
        for (int n = 0; n < n_fpoints; ++n)
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

    // SvdL, 20240901: locations and/or weights of momentum interpolation points are not used (for now), but still passed for consistency with other functions.
    template <typename TF>
    void set_forcing_points_c(
        TF* const restrict tend_c,
        const TF* const restrict fld_c,
        const TF* const restrict fld_u,
        const TF* const restrict fld_v,
        const TF* const restrict fld_w,
        const TF* const restrict boundary_value,
        const int* const restrict gi, const int* const restrict gj, const int* const restrict gk,
        const TF* const restrict rot,
        const int* const restrict ipui, const int* const restrict ipuj, const int* const restrict ipuk, const TF* const restrict c_idw_u,
        const int* const restrict ipvi, const int* const restrict ipvj, const int* const restrict ipvk, const TF* const restrict c_idw_v, 
        const int* const restrict ipwi, const int* const restrict ipwj, const int* const restrict ipwk, const TF* const restrict c_idw_w,
        const int* const restrict ipsi, const int* const restrict ipsj, const int* const restrict ipsk, const TF* const restrict c_idw_s,
        const TF* const restrict db, const TF* const restrict di, const TF* const restrict z0b,
        Boundary_type bc, const TF visc, const int n_fpoints, const int n_idw_loc,
        const int icells, const int ijcells, 
        const double dt)
    {
        const TF  dtf  = static_cast<TF>(dt); // SvdL, 20240901: seems like ugly solution, see how to improve later.

        // For Dirichlet BCs
        if (bc == Boundary_type::Dirichlet_type)
        {
            // Loop over all points to be forced
            for (int n = 0; n < n_fpoints; ++n)
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

    // SvdL, 20240901: locations and/or weights of momentum interpolation points are not used (for now), but still passed for consistency with other functions.
    template <typename TF>
    void set_forcing_points_evisc(
        TF* const restrict fld_evisc,
        const TF* const restrict fld_u,
        const TF* const restrict fld_v,
        const TF* const restrict fld_w,
        // const TF* const restrict boundary_value, //<< SvdL, 20240909: for now not needed for eddy viscosity 
        const int* const restrict gi, const int* const restrict gj, const int* const restrict gk,
        const TF* const restrict rot,
        const int* const restrict ipui, const int* const restrict ipuj, const int* const restrict ipuk, const TF* const restrict c_idw_u,
        const int* const restrict ipvi, const int* const restrict ipvj, const int* const restrict ipvk, const TF* const restrict c_idw_v, 
        const int* const restrict ipwi, const int* const restrict ipwj, const int* const restrict ipwk, const TF* const restrict c_idw_w,
        const int* const restrict ipsi, const int* const restrict ipsj, const int* const restrict ipsk, const TF* const restrict c_idw_s,
        const TF* const restrict db, const TF* const restrict di, const TF* const restrict z0b,
        Boundary_type bc, const TF visc, const int n_fpoints, const int n_idw_loc,
        const int icells, const int ijcells)
    {

        // SvdL, 20240901: IMPLEMENTATION NOTES
        // At immersed boundary, eddy diffusivity is "theoretically" purely driven by the wall-closure model, as we use this same model to set momentum and scalars.
        // Therefore, set K-values accordingly. This implementation requires interpolated momentum at cell center.
        // DO NOT use blending with/extrapolation from LES K-values at second layer (as opposed to Roman et al. [2009] or DeLeon et al. [2018])
        // Blending/extrapolation would make K-value inconsistent with forced momentum/scalars. "Blending" occurs at the cell face in subsequent integration step.
        // This does create a sort-of unsmooth transition on the cell face above.

        // Future: extent with Van Driest (1956) correction and/or make suitable for DNS (see review Verzicco for overview of possibilities, introduces iterative system)
        // DO NOT add molecular viscosity: this is done in the normal diffusion functions.

        // SvdL, 20240918: there is another peculiarity here >> that still needs to be fixed properly <<
        // In case when the grid center coincides with the wall, u,v-momentum is (likely) forced directly to zero value.
        // This is however not allowed for the eddy viscosity, as this would result in zero flux over the boundary, which is needed
        // because the next cell is a "free" grid cell and will need a set K-value for the flux calculation.
        // The correct MOST-consistent value under this conditions can maybe be found from equation... (STILL DO AND IMPLEMENT)

        const int rdim = 9;

        TF u_ip_la;
        TF v_ip_la;
        TF w_ip_la;
        TF umag_ip_la;
        // TF ustar_ip;

        // Loop over all points to be forced
        for (int n = 0; n < n_fpoints; ++n)
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

    template <typename TF>
    void set_ib_points(
        TF* const restrict tend_var,
        const TF* const restrict var,
        const int* const restrict ijk_ib, 
        const TF val,
        const int n_ib,
        const double dt)
    {
        const TF  dtf  = static_cast<TF>(dt); // SvdL, 20240901: seems like ugly solution, see how to improve later.

        for (int nn=0; nn<n_ib; ++nn)
        {
            tend_var[ijk_ib[nn]] = ( val - var[ijk_ib[nn]] ) / dtf;
        }
    }

    template <typename TF>
    void set_ib_points_evisc(
        TF* const restrict var,
        const int* const restrict ijk_ib, 
        const TF val,
        const int n_ib)
    {
        for (int nn=0; nn<n_ib; ++nn)
        {
            var[ijk_ib[nn]] = val;
        }
    }

}

template <typename TF>
Immersed_boundary<TF>::Immersed_boundary(Master &masterin, Grid<TF> &gridin, Fields<TF> &fieldsin, Input &inputin) : master(masterin), grid(gridin), fields(fieldsin),
                                                                                                                     field3d_io(masterin, gridin), boundary_cyclic(masterin, gridin)
{
    // Read IB switch from namelist, and set internal `sw_ib` switch
    std::string sw_ib_str = inputin.get_item<std::string>("IB", "sw_immersed_boundary", "", "0");

    if (sw_ib_str == "0")
        sw_ib = IB_type::Disabled;
    else if (sw_ib_str == "sdf")
        sw_ib = IB_type::SDF;
    else
    {
        std::string error = "\"" + sw_ib_str + "\" is an illegal value for \"sw_ib\"";
        throw std::runtime_error(error);
    }

    if (sw_ib != IB_type::Disabled)
    {
        // Check for use of second order grid and smagorinsky diffusion
        if (grid.get_spatial_order() != Grid_order::Second)
            throw std::runtime_error("Current immersed boundaries only run with second order grids.");

        // if (diff.get_switch() != Diffusion_type::Diff_smag2)
        //     throw std::runtime_error("Current immersed boundaries only run with smagorinsky diffusion.");
        
        // Set a minimum of 3 ghost cells in the horizontal
        const int ijgc = 3;                          // SvdL, 20240828: increased to 3 for increased "options" for projection points
        grid.set_minimum_ghost_cells(ijgc, ijgc, 0); // then interpolation point will be on first layer of ghost cells, and two layers remain for the interpolation

        // Read additional settings
        n_idw_points     = inputin.get_item<int>("IB", "n_idw_points", "", 8);     // SvdL, 20240909: default of 8 seems reasonable..
        n_idw_points_min = inputin.get_item<int>("IB", "n_idw_points_min", "", 4); // SvdL, 20240918: minimum amount of n_idw_points needed for interpolation (4 also reasonable?)

        // Set available masks
        available_masks.insert(available_masks.end(), {"ib"});

        // Read additional parameters from input
        z0bound    = inputin.get_item<TF>("IB", "z0b", "");

        // SvdL, 20240911: should still work...
        // SvdL, 20240724: no check if same scalars are defined in [fields].. they should be AND the order should match
        // If scalars are present (excl. temperature), read in corresponding boundary type
        if (fields.sp.size() > 0)
        {
            swbotlist = inputin.get_list<std::string>("IB", "swbotlist", "");

            if (swbotlist.size() != fields.sp.size())
                throw std::runtime_error("Number of given boundary types does not equal number of scalars.");
        }
    }

}

template <typename TF>
Immersed_boundary<TF>::~Immersed_boundary()
{
}

#ifndef USECUDA
template <typename TF>
void Immersed_boundary<TF>::exec_viscosity()
{
    if (sw_ib == IB_type::Disabled)
        return;

    auto &gd = grid.get_grid_data();

    set_forcing_points_evisc(
        fields.sd.at("evisc")->fld.data(),
        fields.mp.at("u")->fld.data(),
        fields.mp.at("v")->fld.data(),
        fields.mp.at("w")->fld.data(),
        // fpoints.at("s").sbot.at("s").data(),                                                                                      //<< SvdL, 19-06-2023: for now not needed here,  // value of boundary conditions to be enforced
        fpoints.at("s").i.data(), fpoints.at("s").j.data(), fpoints.at("s").k.data(),                                                // points to be forced
        fpoints.at("s").rot.data(),                                                                                                  // rotational matrix for local surface alignment
        fpoints.at("s").ip_u_i.data(), fpoints.at("s").ip_u_j.data(), fpoints.at("s").ip_u_k.data(), fpoints.at("s").c_idw_u.data(), // locations of the neighbouring u-points + weights
        fpoints.at("s").ip_v_i.data(), fpoints.at("s").ip_v_j.data(), fpoints.at("s").ip_v_k.data(), fpoints.at("s").c_idw_v.data(), // locations of the neighbouring v-points + weights
        fpoints.at("s").ip_w_i.data(), fpoints.at("s").ip_w_j.data(), fpoints.at("s").ip_w_k.data(), fpoints.at("s").c_idw_w.data(), // locations of the neighbouring w-points + weights
        fpoints.at("s").ip_s_i.data(), fpoints.at("s").ip_s_j.data(), fpoints.at("s").ip_s_k.data(), fpoints.at("s").c_idw_s.data(), // locations of the neighbouring s-points + weights
        fpoints.at("s").dist_b.data(),                                                                                               // distance nearest immersed boundary point to forcing point
        fpoints.at("s").dist_i.data(),                                                                                               // distance interpolation point to forcing point
        fpoints.at("s").z0b.data(),                                                                                                  // local roughness lengths of forcing points (all scalars will have same for now..)
        Boundary_type::Dirichlet_type,                                                                                               // should contain Boundary_Type:: for all scalars (make variation between scalars possible?), also unused for evisc
        fields.visc, fpoints.at("s").n_fpoints, this->n_idw_points,
        gd.icells, gd.ijcells);

    // eddy viscosity should just be zero inside objects
    set_ib_points_evisc(
        fields.sd.at("evisc")->fld.data(),
        ibpoints.at("s").ijk.data(), TF(0.), 
        ibpoints.at("s").n_ibpoints);

    // SvdL, 20240901: add check if Deardorff scheme is used; then also apply procedure to eviscs? (assume Pr=1 so close to surface?)

    // Enforce cyclic boundary conditions for updated evisc
    boundary_cyclic.exec(fields.sd.at("evisc")->fld.data());
}

template <typename TF>
void Immersed_boundary<TF>::exec(const double dt)
{
    if (sw_ib == IB_type::Disabled)
        return;

    auto &gd = grid.get_grid_data();

    set_forcing_points_u(
        fields.mt.at("u")->fld.data(),
        fields.mt.at("v")->fld.data(),
        fields.mt.at("w")->fld.data(),
        fields.mp.at("u")->fld.data(),
        fields.mp.at("v")->fld.data(),
        fields.mp.at("w")->fld.data(),
        fpoints.at("u").mbot.data(),                                                                                                 // value of boundary conditions to be enforced
        fpoints.at("u").i.data(), fpoints.at("u").j.data(), fpoints.at("u").k.data(),                                                // points to be forced
        fpoints.at("u").rot.data(),                                                                                                  // rotational matrix for local surface alignment
        fpoints.at("u").ip_u_i.data(), fpoints.at("u").ip_u_j.data(), fpoints.at("u").ip_u_k.data(), fpoints.at("u").c_idw_u.data(), // locations of the neighbouring u-points + weights
        fpoints.at("u").ip_v_i.data(), fpoints.at("u").ip_v_j.data(), fpoints.at("u").ip_v_k.data(), fpoints.at("u").c_idw_v.data(), // locations of the neighbouring v-points + weights
        fpoints.at("u").ip_w_i.data(), fpoints.at("u").ip_w_j.data(), fpoints.at("u").ip_w_k.data(), fpoints.at("u").c_idw_w.data(), // locations of the neighbouring w-points + weights
        fpoints.at("u").ip_s_i.data(), fpoints.at("u").ip_s_j.data(), fpoints.at("u").ip_s_k.data(), fpoints.at("u").c_idw_s.data(), // locations of the neighbouring s-points + weights
        fpoints.at("u").dist_b.data(),                                                                                               // distance nearest immersed boundary point to forcing point
        fpoints.at("u").dist_i.data(),                                                                                               // distance interpolation point to forcing point
        fpoints.at("u").z0b.data(),                                                                                                  // local roughness lengths of forcing points
        Boundary_type::Dirichlet_type,                                                                                               // only allow no-slip conditions for momentum (for now..)
        fields.visc, fpoints.at("u").n_fpoints, this->n_idw_points,
        gd.icells, gd.ijcells,
        dt);
    
    set_ib_points(
        fields.mt.at("u")->fld.data(),
        fields.mp.at("u")->fld.data(),
        ibpoints.at("u").ijk.data(), TF(0.), 
        ibpoints.at("u").n_ibpoints,
        dt);

    set_forcing_points_v(
        fields.mt.at("u")->fld.data(),
        fields.mt.at("v")->fld.data(),
        fields.mt.at("w")->fld.data(),
        fields.mp.at("u")->fld.data(),
        fields.mp.at("v")->fld.data(),
        fields.mp.at("w")->fld.data(),
        fpoints.at("v").mbot.data(),                                                                                                 // value of boundary conditions to be enforced
        fpoints.at("v").i.data(), fpoints.at("v").j.data(), fpoints.at("v").k.data(),                                                // points to be forced
        fpoints.at("v").rot.data(),                                                                                                  // rotational matrix for local surface alignment
        fpoints.at("v").ip_u_i.data(), fpoints.at("u").ip_u_j.data(), fpoints.at("v").ip_u_k.data(), fpoints.at("v").c_idw_u.data(), // locations of the neighbouring u-points + weights
        fpoints.at("v").ip_v_i.data(), fpoints.at("v").ip_v_j.data(), fpoints.at("v").ip_v_k.data(), fpoints.at("v").c_idw_v.data(), // locations of the neighbouring v-points + weights
        fpoints.at("v").ip_w_i.data(), fpoints.at("v").ip_w_j.data(), fpoints.at("v").ip_w_k.data(), fpoints.at("v").c_idw_w.data(), // locations of the neighbouring w-points + weights
        fpoints.at("v").ip_s_i.data(), fpoints.at("v").ip_s_j.data(), fpoints.at("v").ip_s_k.data(), fpoints.at("v").c_idw_s.data(), // locations of the neighbouring s-points + weights
        fpoints.at("v").dist_b.data(),                                                                                               // distance nearest immersed boundary point to forcing point
        fpoints.at("v").dist_i.data(),                                                                                               // distance interpolation point to forcing point
        fpoints.at("v").z0b.data(),                                                                                                  // local roughness lengths of forcing points
        Boundary_type::Dirichlet_type,                                                                                               // only allow no-slip conditions for momentum (for now..)
        fields.visc, fpoints.at("v").n_fpoints, this->n_idw_points,
        gd.icells, gd.ijcells, 
        dt);

    set_ib_points(
        fields.mt.at("v")->fld.data(),
        fields.mp.at("v")->fld.data(),
        ibpoints.at("v").ijk.data(), TF(0.), 
        ibpoints.at("v").n_ibpoints,
        dt);

    set_forcing_points_w(
        fields.mt.at("u")->fld.data(),
        fields.mt.at("v")->fld.data(),
        fields.mt.at("w")->fld.data(),
        fields.mp.at("u")->fld.data(),
        fields.mp.at("v")->fld.data(),
        fields.mp.at("w")->fld.data(),
        fpoints.at("w").mbot.data(),                                                                                                 // value of boundary conditions to be enforced
        fpoints.at("w").i.data(), fpoints.at("w").j.data(), fpoints.at("w").k.data(),                                                // points to be forced
        fpoints.at("w").rot.data(),                                                                                                  // rotational matrix for local surface alignment
        fpoints.at("w").ip_u_i.data(), fpoints.at("w").ip_u_j.data(), fpoints.at("w").ip_u_k.data(), fpoints.at("u").c_idw_u.data(), // locations of the neighbouring u-points + weights
        fpoints.at("w").ip_v_i.data(), fpoints.at("w").ip_v_j.data(), fpoints.at("w").ip_v_k.data(), fpoints.at("u").c_idw_v.data(), // locations of the neighbouring v-points + weights
        fpoints.at("w").ip_w_i.data(), fpoints.at("w").ip_w_j.data(), fpoints.at("w").ip_w_k.data(), fpoints.at("u").c_idw_w.data(), // locations of the neighbouring w-points + weights
        fpoints.at("w").ip_s_i.data(), fpoints.at("w").ip_s_j.data(), fpoints.at("w").ip_s_k.data(), fpoints.at("u").c_idw_s.data(), // locations of the neighbouring s-points + weights
        fpoints.at("w").dist_b.data(),                                                                                               // distance nearest immersed boundary point to forcing point
        fpoints.at("w").dist_i.data(),                                                                                               // distance interpolation point to forcing point
        fpoints.at("w").z0b.data(),                                                                                                  // local roughness lengths of forcing points
        Boundary_type::Dirichlet_type,                                                                                               // only allow no-slip conditions for momentum (for now..)
        fields.visc, fpoints.at("w").n_fpoints, this->n_idw_points,
        gd.icells, gd.ijcells,
        dt);

    set_ib_points(
        fields.mt.at("w")->fld.data(),
        fields.mp.at("w")->fld.data(),
        ibpoints.at("w").ijk.data(), TF(0.), 
        ibpoints.at("w").n_ibpoints,
        dt);

    // not required here as we impose ib via tendencies
    // boundary_cyclic.exec(fields.mt.at("u")->fld.data());
    // boundary_cyclic.exec(fields.mt.at("v")->fld.data());
    // boundary_cyclic.exec(fields.mt.at("w")->fld.data());

    for (auto &it : fields.sp)
    {
        set_forcing_points_c(
            fields.st.at(it.first)->fld.data(),
            it.second->fld.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("w")->fld.data(),
            fpoints.at("s").sbot.at(it.first).data(),                                                                                    // value of boundary conditions to be enforced
            fpoints.at("s").i.data(), fpoints.at("s").j.data(), fpoints.at("s").k.data(),                                                // points to be forced
            fpoints.at("s").rot.data(),                                                                                                  // rotational matrix for local surface alignment (although not used yet for scalars: DO calculate, needed for evisc)
            fpoints.at("s").ip_u_i.data(), fpoints.at("s").ip_u_j.data(), fpoints.at("s").ip_u_k.data(), fpoints.at("s").c_idw_u.data(), // locations of the neighbouring u-points + weights
            fpoints.at("s").ip_v_i.data(), fpoints.at("s").ip_v_j.data(), fpoints.at("s").ip_v_k.data(), fpoints.at("s").c_idw_v.data(), // locations of the neighbouring v-points + weights
            fpoints.at("s").ip_w_i.data(), fpoints.at("s").ip_w_j.data(), fpoints.at("s").ip_w_k.data(), fpoints.at("s").c_idw_w.data(), // locations of the neighbouring w-points + weights
            fpoints.at("s").ip_s_i.data(), fpoints.at("s").ip_s_j.data(), fpoints.at("s").ip_s_k.data(), fpoints.at("s").c_idw_s.data(), // locations of the neighbouring s-points + weights
            fpoints.at("s").dist_b.data(),                                                                                               // distance nearest immersed boundary point to forcing point
            fpoints.at("s").dist_i.data(),                                                                                               // distance interpolation point to forcing point
            fpoints.at("s").z0b.data(),                                                                                                  // local roughness lengths of forcing points (all scalars will have same for now..)
            sbc.at(it.first),                                                                                                            // should contain Boundary_Type:: for all scalars (make variation between scalars possible?)
            fields.visc, fpoints.at("s").n_fpoints, this->n_idw_points,
            gd.icells, gd.ijcells,
            dt);

        set_ib_points(
            fields.st.at(it.first)->fld.data(),
            it.second->fld.data(),
            ibpoints.at("s").ijk.data(), TF(0.), 
            ibpoints.at("s").n_ibpoints,
            dt);

    }

    // SvdL, 20240909: plans for much much later... allowing for IB conditions to be updated
}
#endif

template <typename TF>
void Immersed_boundary<TF>::init(Input &inputin, Cross<TF> &cross)
{
    auto &gd = grid.get_grid_data();

    if (sw_ib == IB_type::Disabled)
        return;

    // // SvdL, 20240731: CHECK LATER!! onderstaande moet nog ANDERS!!! DIT WERKT ALLEEN VOOR USER...
    // // TENTATIVE SvdL, 23-07-2024: fix placement of different scalar fields (incl. temperature later). Currently, not expected to work for thermo
    // // DEFINITELY NOT THE NICEST IMPLEMENTATION... I JUST ASSUME ALSO THAT ORDERING OF SCALAR FIELDS IS EQUAL TO READ IN OF SWBOTLIST
    // if (fields.sp.size() > 0)
    // {   
    //     int n = 0;
    //
    //     for (auto &scalar : fields.sp)
    //     {   
    //
    //         std::string swbot = swbotlist[n];
    //
    //         if (swbot == "flux")
    //             sbcbot = Boundary_type::Flux_type;
    //         else if (swbot == "dirichlet")
    //             sbcbot = Boundary_type::Dirichlet_type;
    //         else if (swbot == "neumann")
    //             sbcbot = Boundary_type::Neumann_type;
    //         else
    //         {
    //             std::string error = "IB sbcbot=" + swbot + " is not a valid choice (options: dirichlet, neumann, flux)";
    //             throw std::runtime_error(error);
    //         }    
    //
    //         sbc.emplace(scalar.first, sbcbot);
    //         fpoints.at("s").sbot.emplace(scalar.first, std::vector<TF>()); 
    //
    //         ++n;
    //     }
    // }
    //
    // SvdL, 20240918: not clear about original intention and how to adapt to new version
    // Check input list of cross variables (crosslist)
    // std::vector<std::string>& crosslist_global = cross.get_crosslist();
    // std::vector<std::string>::iterator it = crosslist_global.begin();
    // while (it != crosslist_global.end())
    // {
    //     const std::string fluxbot_ib_string = "fluxbot_ib";
    //     if (has_ending(*it, fluxbot_ib_string))
    //     {
    //         // Strip the ending.
    //         std::string scalar = *it;
    //         scalar.erase(it->length() - fluxbot_ib_string.length());
    //
    //         // Check if array is exists, else cycle.
    //         if (fields.sp.find(scalar) != fields.sp.end())
    //         {
    //             // Remove variable from global list, put in local list
    //             crosslist.push_back(*it);
    //             crosslist_global.erase(it); // erase() returns iterator of next element..
    //         }
    //         else
    //             ++it;
    //     }
    //     else
    //         ++it;
    // }
}

template <typename TF>
void Immersed_boundary<TF>::create(Input &inputin, Netcdf_handle &input_nc)
{

    if (sw_ib == IB_type::Disabled)
        return;

    // Init the toolbox classes.
    boundary_cyclic.init();

    // Get grid and MPI information
    auto& gd  = grid.get_grid_data();
    auto& mpi = master.get_MPI_data();

    const TF no_offset = 0.;
    int nerror = 0;

    // Allocate temperorary fields
    auto sdfs = fields.get_tmp();
    auto tmp1 = fields.get_tmp();
    auto tmp2 = fields.get_tmp();

    // Set structures for the forcing points and ib points (SDF version)
    fpoints.emplace("u", Forcing_points<TF>());
    fpoints.emplace("v", Forcing_points<TF>());
    fpoints.emplace("w", Forcing_points<TF>());
    fpoints.emplace("s", Forcing_points<TF>()); //<< always initialize one for scalars (eddy diffusivity is at this location)

    ibpoints.emplace("u", IB_points<TF>());
    ibpoints.emplace("v", IB_points<TF>());
    ibpoints.emplace("w", IB_points<TF>());
    ibpoints.emplace("s", IB_points<TF>());

    // Read the signed distance fields for momentum positions (temporarily store in allocated tendency fields)
    for (auto &it : fields.mt)
    {
        char filename[256];
        std::sprintf(filename, "sdf.%s.%07d", it.first.c_str(), 0);
        master.print_message("Loading \"%s\" ... ", filename);
        
        if (field3d_io.load_field3d(
                    it.second->fld.data(),
                    tmp1->fld.data(), tmp2->fld.data(),
                    filename, no_offset,
                    gd.kstart, gd.kend))
        {
            master.print_message("FAILED\n");
            ++nerror;
        }
        else
        {
            master.print_message("OK\n");
        }

        // Communicate ghost cells to other MPI subdomains
        boundary_cyclic.exec(it.second->fld.data());

    }

    // Read the signed distance field for scalar positions (temporarily store in tmp-field). Must exist for evisc calculation.
    char filename[256] = "sdf.s.0000000";
    master.print_message("Loading \"%s\" ... ", filename);

    if (field3d_io.load_field3d(
        sdfs->fld.data(),
        tmp1->fld.data(), tmp2->fld.data(),
        filename, no_offset,
        gd.kstart, gd.kend))
    {
        master.print_message("FAILED\n");
        ++nerror;
    }
    else
    {
        master.print_message("OK\n");
    }

    // Communicate ghost cells to other MPI subdomains
    boundary_cyclic.exec(sdfs->fld.data());

    // Now all SDFs are loaded, forcing points, ib points, interpolation and boundary points can be identified
    process_sdf(
            fpoints.at("u"), ibpoints.at("u"),
            fields.mt.at("u")->fld.data(),  // << this contains the SDF
            gd.xh, gd.y, gd.z,              // << location where u is defined
            gd.dx, gd.dy, gd.dzh, gd.dzhi,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.jcells, gd.ijcells,
            this->n_idw_points);

    process_sdf(
            fpoints.at("v"), ibpoints.at("v"),
            fields.mt.at("v")->fld.data(),  // << this contains the SDF
            gd.x, gd.yh, gd.z,              // << location where v is defined
            gd.dx, gd.dy, gd.dzh, gd.dzhi,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.jcells, gd.ijcells,
            this->n_idw_points);

    process_sdf(
            fpoints.at("w"), ibpoints.at("w"),
            fields.mt.at("w")->fld.data(),  // << this contains the SDF
            gd.x, gd.y, gd.zh,              // << location where w is defined
            gd.dx, gd.dy, gd.dz, gd.dzi,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.jcells, gd.ijcells,
            this->n_idw_points);

    process_sdf(
            fpoints.at("s"), ibpoints.at("s"),
            sdfs->fld.data(),               // << this contains the SDF
            gd.x, gd.y, gd.z,               // << location where evisc is defined
            gd.dx, gd.dy, gd.dzh, gd.dzhi,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.jcells, gd.ijcells,
            this->n_idw_points);

    // SvdL, 20240902: temporary quick fix... should be improved later
    // Set roughness of immersed boundaries to constant for all
    fpoints.at("u").z0b.resize(fpoints.at("u").n_fpoints);
    fpoints.at("v").z0b.resize(fpoints.at("v").n_fpoints);
    fpoints.at("w").z0b.resize(fpoints.at("w").n_fpoints);
    fpoints.at("s").z0b.resize(fpoints.at("s").n_fpoints);

    std::fill(fpoints.at("u").z0b.begin(), fpoints.at("u").z0b.end(), z0bound);
    std::fill(fpoints.at("v").z0b.begin(), fpoints.at("v").z0b.end(), z0bound);
    std::fill(fpoints.at("w").z0b.begin(), fpoints.at("w").z0b.end(), z0bound);
    std::fill(fpoints.at("s").z0b.begin(), fpoints.at("s").z0b.end(), z0bound);

    // SvdL, 20240911: temporary quick fix... should be improved later 
    // --> currently UNUSED and internally Dirichlet value 0 is enforced; could be completely removed from functions
    // Set momentum boundary value for momentum forcing points 
    fpoints.at("u").mbot.resize(fpoints.at("u").n_fpoints);
    fpoints.at("v").mbot.resize(fpoints.at("v").n_fpoints);
    fpoints.at("w").mbot.resize(fpoints.at("w").n_fpoints);

    std::fill(fpoints.at("u").mbot.begin(), fpoints.at("u").mbot.end(), TF(0.));
    std::fill(fpoints.at("v").mbot.begin(), fpoints.at("v").mbot.end(), TF(0.));
    std::fill(fpoints.at("w").mbot.begin(), fpoints.at("w").mbot.end(), TF(0.));

    // SvdL, 20240911: this is not the nicest option, default all scalars for now ZERO Dirichlet BC
    // --> this requires that ordering of scalar settings in ini-file are consistent (e.g., in fields, etc.) 
    // --> is there a check for this?
    if (fields.sp.size() > 0)
    {   
        int n = 0;
    
        for (auto &scalar : fields.sp)
        {   
    
            std::string swbot = swbotlist[n];
    
            if (swbot == "flux")
                sbcbot = Boundary_type::Flux_type;
            else if (swbot == "dirichlet")
                sbcbot = Boundary_type::Dirichlet_type;
            else if (swbot == "neumann")
                sbcbot = Boundary_type::Neumann_type;
            else
            {
                std::string error = "IB sbcbot=" + swbot + " is not a valid choice (options: dirichlet, neumann, flux)";
                throw std::runtime_error(error);
            }    
    
            sbc.emplace(scalar.first, sbcbot);
            fpoints.at("s").sbot.emplace(scalar.first, std::vector<TF>()); 
            fpoints.at("s").sbot.at(scalar.first).resize(fpoints.at("s").n_fpoints);

            std::fill(fpoints.at("s").sbot.at(scalar.first).begin(), fpoints.at("s").sbot.at(scalar.first).end(), TF(0.)); //<< SvdL, 20240911: seems to be working (BvS, please confirm?)

            ++n;
        }
    }

    // Setup the interpolation (determine interpolation points an weights) for momentum and scalars
    for (auto &it : fields.mp)
        {
            setup_interpolation(
                    fpoints.at(it.first).xi, fpoints.at(it.first).yi, fpoints.at(it.first).zi,
                    fpoints.at(it.first).ip_u_i, fpoints.at(it.first).ip_u_j, fpoints.at(it.first).ip_u_k, fpoints.at(it.first).c_idw_u, 
                    fpoints.at(it.first).ip_v_i, fpoints.at(it.first).ip_v_j, fpoints.at(it.first).ip_v_k, fpoints.at(it.first).c_idw_v, 
                    fpoints.at(it.first).ip_w_i, fpoints.at(it.first).ip_w_j, fpoints.at(it.first).ip_w_k, fpoints.at(it.first).c_idw_w, 
                    fpoints.at(it.first).ip_s_i, fpoints.at(it.first).ip_s_j, fpoints.at(it.first).ip_s_k, fpoints.at(it.first).c_idw_s, 
                    fpoints.at(it.first).i, fpoints.at(it.first).j, fpoints.at(it.first).k,
                    fpoints.at(it.first).n_fpoints, this->n_idw_points, this->n_idw_points_min,
                    fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(), fields.mt.at("w")->fld.data(), sdfs->fld.data(),                              //<< these contain the SDFs
                    gd.x, gd.y, gd.z, 
                    gd.xh, gd.yh, gd.zh, 
                    gd.dx, gd.dy, 
                    gd.dz, gd.dzh,
                    gd.istart, gd.jstart, gd.kstart,
                    gd.iend, gd.jend, gd.kend,
                    gd.icells, gd.ijcells);
        }

    setup_interpolation(
            fpoints.at("s").xi, fpoints.at("s").yi, fpoints.at("s").zi,
            fpoints.at("s").ip_u_i, fpoints.at("s").ip_u_j, fpoints.at("s").ip_u_k, fpoints.at("s").c_idw_u, 
            fpoints.at("s").ip_v_i, fpoints.at("s").ip_v_j, fpoints.at("s").ip_v_k, fpoints.at("s").c_idw_v, 
            fpoints.at("s").ip_w_i, fpoints.at("s").ip_w_j, fpoints.at("s").ip_w_k, fpoints.at("s").c_idw_w, 
            fpoints.at("s").ip_s_i, fpoints.at("s").ip_s_j, fpoints.at("s").ip_s_k, fpoints.at("s").c_idw_s, 
            fpoints.at("s").i, fpoints.at("s").j, fpoints.at("s").k,
            fpoints.at("s").n_fpoints, this->n_idw_points, this->n_idw_points_min,
            fields.mt.at("u")->fld.data(), fields.mt.at("v")->fld.data(), fields.mt.at("w")->fld.data(), sdfs->fld.data(),                              //<< these contain the SDFs
            gd.x, gd.y, gd.z, 
            gd.xh, gd.yh, gd.zh, 
            gd.dx, gd.dy, 
            gd.dz, gd.dzh,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend, gd.jend, gd.kend,
            gd.icells, gd.ijcells);
 
    // Reset tendency fields and release tmp field (temporary containers for SDFs)
    fields.reset_tendencies();
    fields.release_tmp(sdfs);
    fields.release_tmp(tmp1);
    fields.release_tmp(tmp2);

    // SvdL, 20240824: output some important notes
    master.print_message("SvdL: current statistics between 1st and 2nd forcing layer are likely meaningless at this stage! (Because IB-values are only set after routines in which statistics are calculated\n");
    master.print_message("SvdL: further, the statistics at the forcing points need to be altered to the IBM forcing, with all other tendencies set to zero here.\n");

}

// SvdL, 20240909: fix later.. for now make it an empty function
template <typename TF>
bool Immersed_boundary<TF>::has_mask(std::string mask_name)
{
    // if (std::find(available_masks.begin(), available_masks.end(), mask_name) != available_masks.end())
    //     return true;
    // else
    //     return false;
    return false;
}

// SvdL, 20240909: fix later.. for now make it an empty function
template <typename TF>
void Immersed_boundary<TF>::get_mask(Stats<TF> &stats, std::string mask_name)
{
    auto &gd = grid.get_grid_data();

    // auto mask  = fields.get_tmp();
    // auto maskh = fields.get_tmp();

    // calc_mask(
    //         mask->fld.data(), maskh->fld.data(), dem.data(), gd.z.data(), gd.zh.data(),
    //         gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
    //         gd.icells, gd.ijcells);

    // stats.set_mask_thres("ib", *mask, *maskh, TF(0.5), Stats_mask_type::Plus);

    // fields.release_tmp(mask );
    // fields.release_tmp(maskh);
}

// SvdL, 20240909: fix later..
template <typename TF>
void Immersed_boundary<TF>::exec_cross(Cross<TF> &cross, unsigned long iotime)
{
    auto &gd = grid.get_grid_data();

    // SvdL, 24-05-2023: for now completely disabled.. not sure what it should do exactly.
    // if (cross.get_switch())
    // {
    //     for (auto& s : crosslist)
    //     {
    //         const std::string fluxbot_ib_string = "fluxbot_ib";
    //         if (has_ending(s, fluxbot_ib_string))
    //         {
    //             // Strip the scalar from the fluxbot_ib
    //             std::string scalar = s;
    //             scalar.erase(s.length() - fluxbot_ib_string.length());

    //             auto tmp = fields.get_tmp();

    //             calc_fluxes(
    //                     tmp->flux_bot.data(), k_dem.data(),
    //                     fields.sp.at(scalar)->fld.data(),
    //                     gd.dx,  gd.dy,  gd.dz.data(),
    //                     gd.dxi, gd.dyi, gd.dzhi.data(),
    //                     fields.sp.at(scalar)->visc,
    //                     gd.istart, gd.iend,
    //                     gd.jstart, gd.jend,
    //                     gd.kstart, gd.kend,
    //                     gd.icells, gd.ijcells);

    //             cross.cross_plane(tmp->flux_bot.data(), scalar+"fluxbot_ib", iotime);

    //             fields.release_tmp(tmp);
    //         }
    //     }
    // }
}

template class Immersed_boundary<double>;
template class Immersed_boundary<float>;
