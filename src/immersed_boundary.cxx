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

    // // SvdL, geen idee waar deze voor nodig is...
    // bool has_ending(const std::string& full_string, const std::string& ending)
    // {
    //     if (full_string.length() >= ending.length())
    //         return (0 == full_string.compare(full_string.length() - ending.length(), ending.length(), ending));
    //     else
    //         return false;
    // };

    // void print_statistics(std::vector<int> &ghost_i, std::string name, Master &master)
    // {
    //     int nghost = ghost_i.size();
    //     master.sum(&nghost, 1);

    //     if (master.get_mpiid() == 0)
    //     {
    //         std::string message = "Found: " + std::to_string(nghost) + " IB ghost cells at the " + name + " location";
    //         master.print_message(message);
    //     }
    // }

    // Note 13-06-2023: er moet hier een loop overheen voor alle xi, yi, zi punten..
    template <typename TF>
    void setup_interpolation(
        const TF* const restrict xi_fp, const TF* const restrict yi_fp, const TF* const restrict zi_fp,             // location of interpolation point
        int* const restrict ipui, int* const restrict ipuj, int* const restrict ipuk, TF* const restrict c_idw_u,   // indices of interpolation locations + weights
        int* const restrict ipvi, int* const restrict ipvj, int* const restrict ipvk, TF* const restrict c_idw_v, 
        int* const restrict ipwi, int* const restrict ipwj, int* const restrict ipwk, TF* const restrict c_idw_w,
        int* const restrict ipsi, int* const restrict ipsj, int* const restrict ipsk, TF* const restrict c_idw_s, 
        const std::vector<int>& ijk_fp_u, const std::vector<int>& ijk_ib_u, // these are vectors that contain the ijk combined indices of the forcing points and ib points per grid position
        const std::vector<int>& ijk_fp_v, const std::vector<int>& ijk_ib_v, 
        const std::vector<int>& ijk_fp_w, const std::vector<int>& ijk_ib_w,  
        const std::vector<int>& ijk_fp_s, const std::vector<int>& ijk_ib_s,  
        const int n_fpoints, const int n_idw_loc_min, const int n_idw_loc,
        const std::vector<TF>& x, const std::vector<TF>& y, const std::vector<TF>& z, 
        const std::vector<TF>& xh, const std::vector<TF>& yh, const std::vector<TF>& zh,
        const TF dx, const TF dy, const std::vector<TF>& dz, 
        const int istart, const int jstart, const int kstart,
        const int iend,   const int jend,   const int kend,
        const int icells, const int ijcells,
        const int mpi_offset_x, const int mpi_offset_y)
    {

        for (int nn=0; nn<n_fpoints; ++nn)
        {
            
            const TF xi = xi_fp[nn];
            const TF yi = yi_fp[nn];
            const TF zi = zi_fp[nn];
 
            // Vectors including all neighbours outside IB and not in collection forcing points 
            std::vector<Neighbour<TF>> u_neighbours; // neighbouring u-momentum grid points of current (xi,yi,zi) point
            std::vector<Neighbour<TF>> v_neighbours; // idem for v
            std::vector<Neighbour<TF>> w_neighbours; // idem for w
            std::vector<Neighbour<TF>> s_neighbours; // idem for scalars

            TF dist_max  = TF(0.); 
            TF c_idw_sum = TF(0.);

            // Calculate "semi-nearest" real grid positions of interpolation point (xi,yi,zi) 
            // and subsequently search the neighbourhood for potential interpolation points

            // SvdL, 29-06-2023: this is likely not the best "first" guess of the indices of the interpolation point. 
            // there may be options to get a better initial guess, check if equal to forcing point and then shift based on normal vector direction.
            // But this requires passing even more vectors and more computations.. easier to just extend search region itself..?
            // --> NOTE: visualize on grid representation, floor operator should give position w.r.t. grid center, round operator w.r.t. do grid face
            // const int in_c = static_cast<int>(std::floor(xi/dx)) + istart; 
            // const int jn_c = static_cast<int>(std::floor(yi/dy)) + jstart; 

            // SvdL, 25-07-2024: FORGOT the MPI offset, which is need to make grid search LOCAL
            const int in_c = static_cast<int>(std::round(xi/dx)) + istart - mpi_offset_x; 
            const int jn_c = static_cast<int>(std::round(yi/dy)) + jstart - mpi_offset_y; 
            int kn_c = kstart;

            const int in_f = static_cast<int>(std::round(xi/dx)) + istart - mpi_offset_x; // face position used for u points
            const int jn_f = static_cast<int>(std::round(yi/dy)) + jstart - mpi_offset_y; // face position used for v points
            int kn_f = kstart;                                             // face position used for w points

            for (int k=kstart; k<kend-1; ++k)
            {
                if ( z[k] >= zi ) 
                {
                    if ((z[k]-zi) < (z[k+1]-zi))
                    {
                        kn_c = k; 
                    }
                    else
                    {
                        kn_c = k+1; 
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
                        kn_f = k; 
                    }
                    else
                    {
                        kn_f = k+1; 
                    }
                    break;
                }
            }

            // Limit vertical stencil near surface (interpolation locations may not be below surface)
            const int dk0_c = std::max(-1, kstart-kn_c);
            const int dk0_f = std::max(-1, kstart-kn_f);

            // 1. do for uloc

            dist_max  = 0.; // reset to zero
            c_idw_sum = 0.; // reset to zero

            // Find neighbouring grid points outside IB 
            for (int dk=dk0_c; dk<2; ++dk)
                for (int dj=-1; dj<2; ++dj)
                    for (int di=-1; di<2; ++di)
                    {
                        const int ijk_test = (in_f+di) + (jn_c+dj)*icells + (kn_c+dk)*ijcells; //<< SvdL, 25-07-2024: ijcells is including the ghostcells..
                        bool ijk_in_fp = (std::find(ijk_fp_u.begin(), ijk_fp_u.end(), ijk_test) != ijk_fp_u.end()); //<< SvdL, 25-07-2024: only tests in fp and ib points on current MPI task, should test on all tasks! 
                        bool ijk_in_ib = (std::find(ijk_ib_u.begin(), ijk_ib_u.end(), ijk_test) != ijk_ib_u.end()); //<< THIS IS THEREFORE WRONG... should automatically be solved with signed distance function (test on positive value)

                        // Check if selected grid point is both not a Forcing Point and outside IB.
                        if ( !ijk_in_fp && !ijk_in_ib ) 
                        {
                            const TF distance = absolute_distance(xi, yi, zi, xh[in_f+di], y[jn_c+dj], z[kn_c+dk]);
                            Neighbour<TF> tmp_neighbour = {in_f+di, jn_c+dj, kn_c+dk, distance};
                            u_neighbours.push_back(tmp_neighbour);
                        }
                    }

            // SvdL, 25-07-2024: used for testing, keep for now...
            // if (u_neighbours.empty()) 
            // { 
            //     std::cout << "in_f " << in_f << std::endl; 
            //     std::cout << "jn_c " << jn_c << std::endl; 
            //     std::cout << "kn_c " << kn_c << std::endl; 
            //     // std::cout << "ijk_test " << ijk_test << std::endl; 
            //     std::cout << "ijk_in_fp " << ijk_in_fp << std::endl; 
            //     std::cout << "ijk_in_ib " << ijk_in_ib << std::endl; 
            //     std::cout << "Vector is empty" << std::endl; 
            // } 

            // Sort them on distance
            std::sort(u_neighbours.begin(), u_neighbours.end(), compare_value<TF>);

            // If smallest distance is zero (to within some precision), this point gets weight 1 and the rest is set to zero weight.
            if (u_neighbours[0].distance < TF(1e-7))
            {

                ipui[nn*n_idw_loc]    = u_neighbours[0].i;
                ipuj[nn*n_idw_loc]    = u_neighbours[0].j;
                ipuk[nn*n_idw_loc]    = u_neighbours[0].k;
                c_idw_u[nn*n_idw_loc] = TF(1.);  

                for (int ii=1; ii<n_idw_loc; ++ii)
                {
                    const int in = ii + nn*n_idw_loc; // index gives position in vectors

                    ipui[in]    = in_f;   // just a "FillValue", make sure it carries weight ZERO
                    ipuj[in]    = jn_c;
                    ipuk[in]    = kn_c;
                    c_idw_u[in] = TF(0.); // set weights here to ZERO
                }
            }
            else
            {
                if (u_neighbours.size() < n_idw_loc_min)
                {
                    std::cout << "ERROR: only found " << u_neighbours.size() << " s interpolation points around (less than minimum) ";
                    std::cout << "xi=" << xi << ", yi=" << yi << ", zi=" << zi << std::endl; 
                    throw 1;
                }
                else if (u_neighbours.size() < n_idw_loc)
                {
                    std::cout << "NOTE: only found " << u_neighbours.size() << " s interpolation points around (less than maximum allowed) ";
                    std::cout << "xi=" << xi << ", yi=" << yi << ", zi=" << zi << std::endl; 
                }

                // maximum distance in valid neighbours is available here.
                dist_max = u_neighbours[u_neighbours.size()-1].distance;
                
                for (int ii=0; ii<n_idw_loc; ++ii)
                {
                    const int in = ii + nn*n_idw_loc; // index represent number of Forcing Points for this variable (gives position in vector)
                    if (ii < u_neighbours.size())
                    {   
                        ipui[in]    = u_neighbours[ii].i;
                        ipuj[in]    = u_neighbours[ii].j;
                        ipuk[in]    = u_neighbours[ii].k;
                        c_idw_u[in] = std::pow((dist_max - u_neighbours[ii].distance) / (dist_max * u_neighbours[ii].distance), TF(0.5)) + TF(1e-9); 
                        // SvdL, 12-06-2023: Question: why add the 1e-9 ? And final point always gets weight zero as well (such that it is a wasted point). 
                        // Shouldn't we then increase n_iwd and/or n_idw_loc_max by 1 to account for this?
                        c_idw_sum  += c_idw_u[in];
                    }
                    else
                    {
                        ipui[in]    = in_f;   // just a "FillValue", make sure it carries weight ZERO
                        ipuj[in]    = jn_c;
                        ipuk[in]    = kn_c;
                        c_idw_u[in] = TF(0.); // set weights here to ZERO
                    }
                }
                
                // normalize all weights
                for (int ii=0; ii<n_idw_loc; ++ii)
                {
                    const int in = ii + nn*n_idw_loc; 
                    c_idw_u[in] /= c_idw_sum; 
                }
            }

            // 2. do for vloc

            dist_max  = 0.; // reset to zero
            c_idw_sum = 0.; // reset to zero

            // Find neighbouring grid points outside IB 
            for (int dk=dk0_c; dk<2; ++dk)
                for (int dj=-1; dj<2; ++dj)
                    for (int di=-1; di<2; ++di)
                    {
                        const int ijk_test = (in_c+di) + (jn_f+dj)*icells + (kn_c+dk)*ijcells;
                        bool ijk_in_fp = (std::find(ijk_fp_v.begin(), ijk_fp_v.end(), ijk_test) != ijk_fp_v.end());
                        bool ijk_in_ib = (std::find(ijk_ib_v.begin(), ijk_ib_v.end(), ijk_test) != ijk_ib_v.end());

                        // Check if selected grid point is both outside IB and not a Forcing Point.
                        if ( !ijk_in_fp && !ijk_in_ib )
                        {
                            const TF distance = absolute_distance(xi, yi, zi, x[in_c+di], yh[jn_f+dj], z[kn_c+dk]);
                            Neighbour<TF> tmp_neighbour = {in_c+di, jn_f+dj, kn_c+dk, distance};
                            v_neighbours.push_back(tmp_neighbour);
                        }
                    }

            // Sort them on distance
            std::sort(v_neighbours.begin(), v_neighbours.end(), compare_value<TF>);

            // If smallest distance is zero (to within some precision), this point gets weight 1 and the rest should be zero.
            if (v_neighbours[0].distance < TF(1e-7))
            {
                ipvi[nn*n_idw_loc]    = v_neighbours[0].i;
                ipvj[nn*n_idw_loc]    = v_neighbours[0].j;
                ipvk[nn*n_idw_loc]    = v_neighbours[0].k;
                c_idw_v[nn*n_idw_loc] = TF(1.);  

                for (int ii=1; ii<n_idw_loc; ++ii)
                {
                    const int in = ii + nn*n_idw_loc; 
                    ipvi[in]    = in_c;   // just a "FillValue", make sure it carries weight ZERO
                    ipvj[in]    = jn_f;
                    ipvk[in]    = kn_c;
                    c_idw_v[in] = TF(0.); // set weights here to ZERO
                }
            }
            else
            {
                if (v_neighbours.size() < n_idw_loc_min)
                {
                    std::cout << "ERROR: only found " << v_neighbours.size() << " s interpolation points around (less than minimum) ";
                    std::cout << "xi=" << xi << ", yi=" << yi << ", zi=" << zi << std::endl; 
                    throw 1;
                }
                else if (v_neighbours.size() < n_idw_loc)
                {
                    std::cout << "NOTE: only found " << v_neighbours.size() << " s interpolation points around (less than maximum allowed) ";
                    std::cout << "xi=" << xi << ", yi=" << yi << ", zi=" << zi << std::endl; 
                }

                // maximum distance in valid neighbours is available here.
                dist_max = v_neighbours[v_neighbours.size()-1].distance;
                
                for (int ii=0; ii<n_idw_loc; ++ii)
                {
                    const int in = ii + nn*n_idw_loc; // index represent number of Forcing Points for this variable (gives position in vector)
                    if (ii < v_neighbours.size())
                    {   
                        ipvi[in]    = v_neighbours[ii].i;
                        ipvj[in]    = v_neighbours[ii].j;
                        ipvk[in]    = v_neighbours[ii].k;
                        c_idw_v[in] = std::pow((dist_max - v_neighbours[ii].distance) / (dist_max * v_neighbours[ii].distance), TF(0.5)) + TF(1e-9); 
                        c_idw_sum  += c_idw_v[in];
                    }
                    else
                    {
                        ipvi[in]    = in_c;   // just a "FillValue", make sure it carries weight ZERO
                        ipvj[in]    = jn_f;
                        ipvk[in]    = kn_c;
                        c_idw_v[in] = TF(0.); // set weights here to ZERO
                    }
                }

                // normalize all weights
                for (int ii=0; ii<n_idw_loc; ++ii)
                {
                    const int in = ii + nn*n_idw_loc; 
                    c_idw_v[in] /= c_idw_sum; 
                }
            }

            // 3. do for wloc

            dist_max  = 0.; // reset to zero
            c_idw_sum = 0.; // reset to zero

            // Find neighbouring grid points outside IB 
            for (int dk=dk0_f; dk<2; ++dk)
                for (int dj=-1; dj<2; ++dj)
                    for (int di=-1; di<2; ++di)
                    {
                        const int ijk_test = (in_c+di) + (jn_c+dj)*icells + (kn_f+dk)*ijcells;
                        bool ijk_in_fp = (std::find(ijk_fp_w.begin(), ijk_fp_w.end(), ijk_test) != ijk_fp_w.end()); 
                        bool ijk_in_ib = (std::find(ijk_ib_w.begin(), ijk_ib_w.end(), ijk_test) != ijk_ib_w.end());

                        // Check if selected grid point is both outside IB and not a Forcing Point.
                        if ( !ijk_in_fp && !ijk_in_ib )
                        {
                            const TF distance = absolute_distance(xi, yi, zi, x[in_c+di], y[jn_c+dj], zh[kn_f+dk]);
                            Neighbour<TF> tmp_neighbour = {in_c+di, jn_c+dj, kn_f+dk, distance};
                            w_neighbours.push_back(tmp_neighbour);
                        }
                    }

            // Sort them on distance
            std::sort(w_neighbours.begin(), w_neighbours.end(), compare_value<TF>);
            
            // If smallest distance is zero (to within some precision), this point gets weight 1 and the rest should be zero.
            if (w_neighbours[0].distance < TF(1e-7))
            {
                ipwi[nn*n_idw_loc]    = w_neighbours[0].i;
                ipwj[nn*n_idw_loc]    = w_neighbours[0].j;
                ipwk[nn*n_idw_loc]    = w_neighbours[0].k;
                c_idw_w[nn*n_idw_loc] = TF(1.);  

                for (int ii=1; ii<n_idw_loc; ++ii)
                {
                    const int in = ii + nn*n_idw_loc;

                    ipwi[in]    = in_c;   // just a "FillValue", make sure it carries weight ZERO
                    ipwj[in]    = jn_c;
                    ipwk[in]    = kn_f;
                    c_idw_w[in] = TF(0.); // set weights here to ZERO
                }
            }
            else
            {
                if (w_neighbours.size() < n_idw_loc_min)
                {
                    std::cout << "ERROR: only found " << w_neighbours.size() << " s interpolation points around (less than minimum) ";
                    std::cout << "xi=" << xi << ", yi=" << yi << ", zi=" << zi << std::endl; 
                    throw 1;
                }
                else if (w_neighbours.size() < n_idw_loc)
                {
                    std::cout << "NOTE: only found " << w_neighbours.size() << " s interpolation points around (less than maximum allowed) ";
                    std::cout << "xi=" << xi << ", yi=" << yi << ", zi=" << zi << std::endl; 
                }

                // maximum distance in valid neighbours is available here.
                dist_max = w_neighbours[w_neighbours.size()-1].distance;  
                
                for (int ii=0; ii<n_idw_loc; ++ii)
                {
                    const int in = ii + nn*n_idw_loc;
                    if (ii < w_neighbours.size())
                    {   
                        ipwi[in]    = w_neighbours[ii].i;
                        ipwj[in]    = w_neighbours[ii].j;
                        ipwk[in]    = w_neighbours[ii].k;
                        c_idw_w[in] = std::pow((dist_max - w_neighbours[ii].distance) / (dist_max * w_neighbours[ii].distance), TF(0.5)) + TF(1e-9); 
                        c_idw_sum  += c_idw_w[in];
                    }
                    else
                    {
                        ipwi[in]    = in_c;   // just a "FillValue", make sure it carries weight ZERO
                        ipwj[in]    = jn_c;
                        ipwk[in]    = kn_f;
                        c_idw_w[in] = TF(0.); // set weights here to ZERO
                    }
                }

                // normalize all weights
                for (int ii=0; ii<n_idw_loc; ++ii)
                {
                    const int in = ii + nn*n_idw_loc; 
                    c_idw_w[in] /= c_idw_sum; 
                }
            }

            // 4. do for sloc

            dist_max  = 0.; // reset to zero
            c_idw_sum = 0.; // reset to zero

            // Find neighbouring grid points outside IB 
            for (int dk=dk0_c; dk<2; ++dk)
                for (int dj=-1; dj<2; ++dj)
                    for (int di=-1; di<2; ++di)
                    {
                        const int ijk_test = (in_c+di) + (jn_c+dj)*icells + (kn_c+dk)*ijcells;
                        bool ijk_in_fp = (std::find(ijk_fp_s.begin(), ijk_fp_s.end(), ijk_test) != ijk_fp_s.end()); 
                        bool ijk_in_ib = (std::find(ijk_ib_s.begin(), ijk_ib_s.end(), ijk_test) != ijk_ib_s.end());

                        // Check if selected grid point is both outside IB and not a Forcing Point.
                        if ( !ijk_in_fp && !ijk_in_ib ) 
                        {
                            const TF distance = absolute_distance(xi, yi, zi, x[in_c+di], y[jn_c+dj], z[kn_c+dk]);
                            Neighbour<TF> tmp_neighbour = {in_c+di, jn_c+dj, kn_c+dk, distance};
                            s_neighbours.push_back(tmp_neighbour);
                        }
                    }

            // Sort them on distance
            std::sort(s_neighbours.begin(), s_neighbours.end(), compare_value<TF>);
            
            // If smallest distance is zero (to within some precision), this point gets weight 1 and the rest should be zero.
            if (s_neighbours[0].distance < TF(1e-7))
            {
                ipsi[nn*n_idw_loc]    = s_neighbours[0].i;
                ipsj[nn*n_idw_loc]    = s_neighbours[0].j;
                ipsk[nn*n_idw_loc]    = s_neighbours[0].k;
                c_idw_s[nn*n_idw_loc] = TF(1.);  

                for (int ii=1; ii<n_idw_loc; ++ii)
                {
                    const int in = ii + nn*n_idw_loc; 

                    ipsi[in]    = in_c;   // just a "FillValue", make sure it carries weight ZERO
                    ipsj[in]    = jn_c;
                    ipsk[in]    = kn_c;
                    c_idw_s[in] = TF(0.); // set weights here to ZERO
                }
            }
            else
            {
                if (s_neighbours.size() < n_idw_loc_min) 
                {
                    std::cout << "ERROR: only found " << s_neighbours.size() << " s interpolation points around (less than minimum) ";
                    std::cout << "xi=" << xi << ", yi=" << yi << ", zi=" << zi << std::endl; 
                    throw 1;
                }
                else if (s_neighbours.size() < n_idw_loc)
                {
                    std::cout << "NOTE: only found " << s_neighbours.size() << " s interpolation points around (less than maximum allowed) ";
                    std::cout << "xi=" << xi << ", yi=" << yi << ", zi=" << zi << std::endl; 
                }

                // maximum distance in valid neighbours is available here.
                dist_max = s_neighbours[s_neighbours.size()-1].distance;  
                
                for (int ii=0; ii<n_idw_loc; ++ii)
                {
                    const int in = ii + nn*n_idw_loc; 
                    if (ii <s_neighbours.size())
                    {   
                        ipsi[in]    = s_neighbours[ii].i;
                        ipsj[in]    = s_neighbours[ii].j;
                        ipsk[in]    = s_neighbours[ii].k;
                        c_idw_s[in] = std::pow((dist_max - s_neighbours[ii].distance) / (dist_max * s_neighbours[ii].distance), TF(0.5)) + TF(1e-9); 
                        c_idw_sum  += c_idw_s[in];
                    }
                    else
                    {
                        ipsi[in]    = in_c;   // just a "FillValue", make sure it carries weight ZERO
                        ipsj[in]    = jn_c;
                        ipsk[in]    = kn_c;
                        c_idw_s[in] = TF(0.); // set weights here to ZERO
                    }
                }

                // normalize all weights
                for (int ii=0; ii<n_idw_loc; ++ii)
                {
                    const int in = ii + nn*n_idw_loc; 
                    c_idw_s[in] /= c_idw_sum; 
                }
            }
        }
    }   

    template <typename TF>
    void calculate_rotation_matrix(
        TF* const restrict rot,
        const TF* const restrict norm, 
        const int n_fpoints)
    {
        const int rdim = 9;

         // Loop over all points to be forced
        for (int n = 0; n < n_fpoints; ++n)
        {

            // seperate out the x,y,z components of the normal vector.
            const TF normx = norm[3*n];
            const TF normy = norm[3*n+1];
            const TF normz = norm[3*n+2];

            // SvdL, 29-06-2023: still add description of how this matrix is set up.
            // writing the whole system out likely makes rotation matrix redundant at all, as operations can easily be done inside the functions.
            // also check for further degenerate cases, matrix below could maybe go wrong with precision or float to double conversions.
    
            if (normx == TF(1.) && normy == TF(0.) && normz == TF(0.))
            {
                rot[rdim * n]     = TF(0.);                              
                rot[rdim * n + 1] = TF(0.);
                rot[rdim * n + 2] = TF(-1.);
                rot[rdim * n + 3] = TF(0.);
                rot[rdim * n + 4] = TF(1.);
                rot[rdim * n + 5] = TF(0.);
                rot[rdim * n + 6] = TF(1.);
                rot[rdim * n + 7] = TF(0.);
                rot[rdim * n + 8] = TF(0.);
            }
            else if (normx == TF(-1.) && normy == TF(0.) && normz == TF(0.))
            {
                rot[rdim * n]     = TF(0.);                              
                rot[rdim * n + 1] = TF(0.);
                rot[rdim * n + 2] = TF(1.);
                rot[rdim * n + 3] = TF(0.);
                rot[rdim * n + 4] = TF(1.);
                rot[rdim * n + 5] = TF(0.);
                rot[rdim * n + 6] = TF(-1.);
                rot[rdim * n + 7] = TF(0.);
                rot[rdim * n + 8] = TF(0.);
            }
            else
            {
                const TF scale = std::sqrt( fm::pow2(normy) + fm::pow2(normz)); 

                rot[rdim * n]     = (fm::pow2(normy) + fm::pow2(normz)) / scale;                              
                rot[rdim * n + 1] = -normx*normy/scale;
                rot[rdim * n + 2] = -normx*normz/scale;
                rot[rdim * n + 3] = TF(0.);
                rot[rdim * n + 4] = normz/scale;
                rot[rdim * n + 5] = -normy/scale;
                rot[rdim * n + 6] = normx;
                rot[rdim * n + 7] = normy;
                rot[rdim * n + 8] = normz;
            }
        }

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
        const int* const restrict ipsi, const int* const restrict ipsj, const int* const restrict ipsk, const TF* const restrict c_idw_s, // SvdL, 08-06-2023: not used for now..
        const TF* const restrict db, const TF* const restrict di, const TF* const restrict z0b,
        Boundary_type bc, const TF visc, const int n_fpoints, const int n_idw_loc,
        const int icells, const int ijcells,
        const double dt)
    {

        const int rdim = 9;                                        
        const TF  dtf  = static_cast<TF>(dt);// SvdL, 15-06-2023: just an intermediate solution/test for now. see how to fix better later.

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
                const int ii = i + n * n_idw_loc;                                   // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE
                const int ijku = ipui[ii] + ipuj[ii] * icells + ipuk[ii] * ijcells; // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE
                const int ijkv = ipvi[ii] + ipvj[ii] * icells + ipvk[ii] * ijcells; // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE
                const int ijkw = ipwi[ii] + ipwj[ii] * icells + ipwk[ii] * ijcells; // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE
                
                // Do the correction based on the auxiliary velocity (i.e. intermediate velocity at next timestep without pressure forcing).
                u_ip += c_idw_u[ii] * (fld_u[ijku] + dtf * tend_u[ijku]);
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
                u_fp_la = u_ip_la * std::log(db[n] / z0b[n]) / std::log((db[n] + di[n]) / z0b[n]);
                v_fp_la = v_ip_la * std::log(db[n] / z0b[n]) / std::log((db[n] + di[n]) / z0b[n]);
                w_fp_la = w_ip_la * fm::pow2(db[n] / (db[n] + di[n]));

                // 4. rotate back to standard grid alginment (only one component is needed here), 
                // AND overwrite old tendency at forcing point with new one to achieve this.
                tend_u[ijkf] = ( (r11 * u_fp_la + r21 * v_fp_la + r31 * w_ip_la) - fld_u[ijkf] ) / dtf;
            }
            else // SvdL, 29-06-2023: change later into Van Driest like correction..
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
        const int* const restrict ipsi, const int* const restrict ipsj, const int* const restrict ipsk, const TF* const restrict c_idw_s, // SvdL, 08-06-2023: not used for now..
        const TF* const restrict db, const TF* const restrict di, const TF* const restrict z0b,
        Boundary_type bc, const TF visc, const int n_fpoints, const int n_idw_loc,
        const int icells, const int ijcells,
        const double dt)
    {

        const int rdim = 9;                                       
        const TF  dtf  = static_cast<TF>(dt); // SvdL, 15-06-2023: just an intermediate solution/test for now. see how to fix better later.

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
                const int ii = i + n * n_idw_loc;                                   // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE
                const int ijku = ipui[ii] + ipuj[ii] * icells + ipuk[ii] * ijcells; // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE
                const int ijkv = ipvi[ii] + ipvj[ii] * icells + ipvk[ii] * ijcells; // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE
                const int ijkw = ipwi[ii] + ipwj[ii] * icells + ipwk[ii] * ijcells; // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE

                // Do the correction based on the auxiliary velocity (i.e. intermediate velocity at next timestep without pressure forcing).
                u_ip += c_idw_u[ii] * (fld_u[ijku] + dtf * tend_u[ijku]);
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
                u_fp_la = u_ip_la * std::log(db[n] / z0b[n]) / std::log((db[n] + di[n]) / z0b[n]);
                v_fp_la = v_ip_la * std::log(db[n] / z0b[n]) / std::log((db[n] + di[n]) / z0b[n]);
                w_fp_la = w_ip_la * fm::pow2(db[n] / (db[n] + di[n]));

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
        const int* const restrict ipsi, const int* const restrict ipsj, const int* const restrict ipsk, const TF* const restrict c_idw_s, // SvdL, 08-06-2023: not used for now..
        const TF* const restrict db, const TF* const restrict di, const TF* const restrict z0b,
        Boundary_type bc, const TF visc, const int n_fpoints, const int n_idw_loc,
        const int icells, const int ijcells,
        const double dt)
    {

        const int rdim = 9;                                        
        const TF  dtf  = static_cast<TF>(dt); // SvdL, 15-06-2023: just an intermediate solution/test for now. see how to fix better later.

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
                const int ii = i + n * n_idw_loc;                                   // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE
                const int ijku = ipui[ii] + ipuj[ii] * icells + ipuk[ii] * ijcells; // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE
                const int ijkv = ipvi[ii] + ipvj[ii] * icells + ipvk[ii] * ijcells; // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE
                const int ijkw = ipwi[ii] + ipwj[ii] * icells + ipwk[ii] * ijcells; // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE

                // Do the correction based on the auxiliary velocity (i.e. intermediate velocity at next timestep without pressure forcing).
                u_ip += c_idw_u[ii] * (fld_u[ijku] + dtf * tend_u[ijku]);
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
                u_fp_la = u_ip_la * std::log(db[n] / z0b[n]) / std::log((db[n] + di[n]) / z0b[n]);
                v_fp_la = v_ip_la * std::log(db[n] / z0b[n]) / std::log((db[n] + di[n]) / z0b[n]);
                w_fp_la = w_ip_la * fm::pow2(db[n] / (db[n] + di[n]));

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

    // SvdL, 12-06-2023: locations and or weights of momentum interpolation points are not used (for now), but still passed for consistency with other functions.
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
        const TF  dtf  = static_cast<TF>(dt); // SvdL, 15-06-2023: just an intermediate solution/test for now. see how to fix better later.

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
                    const int ii = i + n * n_idw_loc;                                   // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE
                    const int ijki = ipsi[ii] + ipsj[ii] * icells + ipsk[ii] * ijcells; // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE
                    c_ip += c_idw_s[ii] * (fld_c[ijki] + dtf * tend_c[ijki] );
                }

                // 2. calculate scalar at forcing point and force at once
                // for now, (1) neglect stability effects (requires "fine enough" grid), and (2) assume both points are in logarithmic layer
                if (db[n] > z0b[n])
                {
                    tend_c[ijkf] = ( ( (c_ip - boundary_value[n]) * std::log(db[n] / z0b[n]) / std::log((db[n] + di[n]) / z0b[n]) + boundary_value[n] ) - fld_c[ijkf] ) / dtf;
                }
                else
                {
                    tend_c[ijkf] = ( boundary_value[n] - fld_c[ijkf] ) / dtf;
                }
            }
        }
        else if (bc == Boundary_type::Flux_type)
        {
            return; // SvdL, 29-06-2023: still implement later
        }
        else if (bc == Boundary_type::Neumann_type)
        {
            return; // SvdL, 29-06-2023: still implement later
        }

    }

    // SvdL, 12-06-2023: locations and or weights of momentum interpolation points are not used (for now), but still passed for consistency with other functions.
    template <typename TF>
    void set_forcing_points_evisc(
        TF* const restrict fld_evisc,
        const TF* const restrict fld_u,
        const TF* const restrict fld_v,
        const TF* const restrict fld_w,
        // const TF* const restrict boundary_value, //<< SvdL, 19-06-2023: for now not needed for eddy viscosity 
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

        // SvdL, 21-05-2023: IMPLEMENTATION NOTES
        // At immersed boundary, eddy diffusivity is "theoretically" purely driven by the wall-closure model, as we use this same model to set momentum and scalars.
        // Therefore, set K-values accordingly. This implementation requires interpolated momentum at cell center.
        // DO NOT use blending with/extrapolation from LES K-values at second layer (as opposed to Roman et al. [2009] or DeLeon et al. [2018])
        // Blending/extrapolation would make K-value inconsistent with forced momentum/scalars. "Blending" occurs at the cell face in subsequent integration step.
        // This does create a sort-of unsmooth transition on the cell face above.

        // Future: extent with Van Driest (1956) correction and/or make suitable for DNS.
        // DO NOT add molecular viscosity: this is done in the normal diffusion functions.

        // SvdL, 23-05-2023: there is another peculiarity here >> that still needs to be fixed properly <<
        // In case when the grid center coincides with the wall, u,v-momentum is (likely) forced directly to zero value.
        // This is however not allowed for the eddy viscosity, as this would result in zero flux over the boundary, which is needed
        // because the next cell is a "free" grid cell and will need a set K-value for the flux calculation.
        // The correct MOST-consistent value under this conditions can be found from equation... (STILL DO AND IMPLEMENT)

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
                const int ii = i + n * n_idw_loc;                                   // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE
                const int ijku = ipui[ii] + ipuj[ii] * icells + ipuk[ii] * ijcells; // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE
                const int ijkv = ipvi[ii] + ipvj[ii] * icells + ipvk[ii] * ijcells; // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE
                const int ijkw = ipwi[ii] + ipwj[ii] * icells + ipwk[ii] * ijcells; // SvdL, 20-05-2023: CHECK LOCATIONS!! NOTE

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
                fld_evisc[ijkf] = fm::pow2(Constants::kappa<TF>) * umag_ip_la * db[n] / std::log((db[n] + di[n]) / z0b[n]);
            }
            else
            {
                // SvdL, 24-05-2023: Currently defaulting to WRONG evisc-value due to zero boundary distance..
                fld_evisc[ijkf] = fm::pow2(Constants::kappa<TF>) * umag_ip_la * z0b[n] / std::log((db[n] + di[n]) / z0b[n]);
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
        const TF  dtf  = static_cast<TF>(dt);// SvdL, 15-06-2023: just an intermediate solution/test for now. see how to fix better later.

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
    // else if (sw_ib_str == "obj") // SvdL, 21-05-2023: potentially implement others options later, e.g., OBJ/STL/SDF 
    //     sw_ib = IB_type::OBJ;    // (practically the same as by user defined input but now all processing should be moved from python to microhh)
    else if (sw_ib_str == "user")   // SvdL, in User defined input: specificy e.g. via netcdf the fpoints-data directly (from preprocess in python)
        sw_ib = IB_type::User;
    else
    {
        std::string error = "\"" + sw_ib_str + "\" is an illegal value for \"sw_ib\"";
        throw std::runtime_error(error);
    }

    // SvdL, 20-05-2023: put check on second order grids and use of smagorinsky diffusion
    // SvdL, working smagorinsky should easily be extendable to DNS in 2nd order...
    if (grid.get_spatial_order() != Grid_order::Second)
        throw std::runtime_error("Current immersed boundaries only run with second order grids.");

    // if (diff.get_switch() != Diffusion_type::Diff_smag2)
    //     throw std::runtime_error("Current immersed boundaries only run with smagorinsky diffusion.");

    if (sw_ib != IB_type::Disabled)
    {
        // Set a minimum of 3 ghost cells in the horizontal
        const int ijgc = 3;                          // SvdL, 01-06-2023: increased this to three to cover for possibility that fpoint is on MPI boundary, 
        grid.set_minimum_ghost_cells(ijgc, ijgc, 0); // then interpolation point will be on first layer of ghost cells, and two layers remain for the interpolation

        // Read additional settings
        n_idw_points     = inputin.get_item<int>("IB", "n_idw_points", "", 8);     // SvdL, 21-05-2023: default of 8 seems reasonable.. saw this in paper Lundquist et al.
        n_idw_points_min = inputin.get_item<int>("IB", "n_idw_points_min", "", 4); // SvdL, 06-06-2023: minimum amount of n_idw_points needed for interpolation (NOTE, still check amount)

        // Set available masks
        available_masks.insert(available_masks.end(), {"ib"});

        // SvdL, 24-07-2024: no check if same scalars are defined in [fields].. they should be AND the order should match
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
        // fpoints.at("s").sbot.at("s").data(),                                                         //<< SvdL, 19-06-2023: for now not needed here,  // value of boundary conditions to be enforced
        fpoints.at("s").i.data(), fpoints.at("s").j.data(), fpoints.at("s").k.data(),                   // points to be forced
        fpoints.at("s").rot.data(),                                                                     // rotational matrix for local surface alignment
        fpoints.at("s").ip_u_i.data(), fpoints.at("s").ip_u_j.data(), fpoints.at("s").ip_u_k.data(), fpoints.at("s").c_idw_u.data(), // locations of the neighbouring u-points + weights
        fpoints.at("s").ip_v_i.data(), fpoints.at("s").ip_v_j.data(), fpoints.at("s").ip_v_k.data(), fpoints.at("s").c_idw_v.data(), // locations of the neighbouring v-points + weights
        fpoints.at("s").ip_w_i.data(), fpoints.at("s").ip_w_j.data(), fpoints.at("s").ip_w_k.data(), fpoints.at("s").c_idw_w.data(), // locations of the neighbouring w-points + weights
        fpoints.at("s").ip_s_i.data(), fpoints.at("s").ip_s_j.data(), fpoints.at("s").ip_s_k.data(), fpoints.at("s").c_idw_s.data(), // locations of the neighbouring s-points + weights
        fpoints.at("s").db.data(),                                                                      // distance nearest immersed boundary point to forcing point
        fpoints.at("s").di.data(),                                                                      // distance interpolation point to forcing point
        fpoints.at("s").z0b.data(),                                                                     // local roughness lengths of forcing points (all scalars will have same for now..)
        Boundary_type::Dirichlet_type,                                                                  // should contain Boundary_Type:: for all scalars (make variation between scalars possible?), also unused for evisc
        fields.visc, fpoints.at("s").i.size(), n_idw_points,
        gd.icells, gd.ijcells);

    // eddy viscosity should just be zero inside..
    set_ib_points_evisc(
        fields.sd.at("evisc")->fld.data(),
        fpoints.at("s").ijk_ib.data(), TF(0.), 
        fpoints.at("s").ijk_ib.size());

    // SvdL, 25-06-2023: STILL ADD FOR EVISCS, if needed?? -> Deardorff scheme?

    // this one is needed.. SvdL, 15-06-2023: do check 
    boundary_cyclic.exec(fields.sd.at("evisc")->fld.data());
}

template <typename TF>
void Immersed_boundary<TF>::exec(const double dt)
{
    if (sw_ib == IB_type::Disabled)
        return;

    auto &gd = grid.get_grid_data();

    // SvdL, 21-05-2023: this approach is only valid IF forcing point of one variable does not end up in interpolation component of another.
    // NOTE: 14-06-2023: waarom mag je hier ip_u_i data gewoon als array meegeven, (zie function definition) maar moet je hem ergens anders bijv weer als vector behandelen..
    set_forcing_points_u(
        fields.mt.at("u")->fld.data(),
        fields.mt.at("v")->fld.data(),
        fields.mt.at("w")->fld.data(),
        fields.mp.at("u")->fld.data(),
        fields.mp.at("v")->fld.data(),
        fields.mp.at("w")->fld.data(),
        fpoints.at("u").mbot.data(),                                                                                            // value of boundary conditions to be enforced
        fpoints.at("u").i.data(), fpoints.at("u").j.data(), fpoints.at("u").k.data(),                                           // points to be forced
        fpoints.at("u").rot.data(),                                                                                             // rotational matrix for local surface alignment
        fpoints.at("u").ip_u_i.data(), fpoints.at("u").ip_u_j.data(), fpoints.at("u").ip_u_k.data(), fpoints.at("u").c_idw_u.data(), // locations of the neighbouring u-points + weights
        fpoints.at("u").ip_v_i.data(), fpoints.at("u").ip_v_j.data(), fpoints.at("u").ip_v_k.data(), fpoints.at("u").c_idw_v.data(), // locations of the neighbouring v-points + weights
        fpoints.at("u").ip_w_i.data(), fpoints.at("u").ip_w_j.data(), fpoints.at("u").ip_w_k.data(), fpoints.at("u").c_idw_w.data(), // locations of the neighbouring w-points + weights
        fpoints.at("u").ip_s_i.data(), fpoints.at("u").ip_s_j.data(), fpoints.at("u").ip_s_k.data(), fpoints.at("u").c_idw_s.data(), // locations of the neighbouring s-points + weights
        fpoints.at("u").db.data(),                                                                                              // distance nearest immersed boundary point to forcing point
        fpoints.at("u").di.data(),                                                                                              // distance interpolation point to forcing point
        fpoints.at("u").z0b.data(),                                                                                             // local roughness lengths of forcing points
        Boundary_type::Dirichlet_type,                                                                                          // only allow no-slip conditions for momentum (for now..)
        fields.visc, fpoints.at("u").i.size(), n_idw_points,
        gd.icells, gd.ijcells,
        dt);
    
    set_ib_points(
        fields.mt.at("u")->fld.data(),
        fields.mp.at("u")->fld.data(),
        fpoints.at("u").ijk_ib.data(), TF(0.), 
        fpoints.at("u").ijk_ib.size(),
        dt);

    set_forcing_points_v(
        fields.mt.at("u")->fld.data(),
        fields.mt.at("v")->fld.data(),
        fields.mt.at("w")->fld.data(),
        fields.mp.at("u")->fld.data(),
        fields.mp.at("v")->fld.data(),
        fields.mp.at("w")->fld.data(),
        fpoints.at("v").mbot.data(),                                                                    // value of boundary conditions to be enforced
        fpoints.at("v").i.data(), fpoints.at("v").j.data(), fpoints.at("v").k.data(),                   // points to be forced
        fpoints.at("v").rot.data(),                                                                     // rotational matrix for local surface alignment
        fpoints.at("v").ip_u_i.data(), fpoints.at("u").ip_u_j.data(), fpoints.at("v").ip_u_k.data(), fpoints.at("v").c_idw_u.data(), // locations of the neighbouring u-points + weights
        fpoints.at("v").ip_v_i.data(), fpoints.at("v").ip_v_j.data(), fpoints.at("v").ip_v_k.data(), fpoints.at("v").c_idw_v.data(), // locations of the neighbouring v-points + weights
        fpoints.at("v").ip_w_i.data(), fpoints.at("v").ip_w_j.data(), fpoints.at("v").ip_w_k.data(), fpoints.at("v").c_idw_w.data(), // locations of the neighbouring w-points + weights
        fpoints.at("v").ip_s_i.data(), fpoints.at("v").ip_s_j.data(), fpoints.at("v").ip_s_k.data(), fpoints.at("v").c_idw_s.data(), // locations of the neighbouring s-points + weights
        fpoints.at("v").db.data(),                                                                      // distance nearest immersed boundary point to forcing point
        fpoints.at("v").di.data(),                                                                      // distance interpolation point to forcing point
        fpoints.at("v").z0b.data(),                                                                     // local roughness lengths of forcing points
        Boundary_type::Dirichlet_type,                                                                  // only allow no-slip conditions for momentum (for now..)
        fields.visc, fpoints.at("v").i.size(), n_idw_points,
        gd.icells, gd.ijcells, 
        dt);

    set_ib_points(
        fields.mt.at("v")->fld.data(),
        fields.mp.at("v")->fld.data(),
        fpoints.at("v").ijk_ib.data(), TF(0.), 
        fpoints.at("v").ijk_ib.size(),
        dt);

    set_forcing_points_w(
        fields.mt.at("u")->fld.data(),
        fields.mt.at("v")->fld.data(),
        fields.mt.at("w")->fld.data(),
        fields.mp.at("u")->fld.data(),
        fields.mp.at("v")->fld.data(),
        fields.mp.at("w")->fld.data(),
        fpoints.at("w").mbot.data(),                                                                    // value of boundary conditions to be enforced
        fpoints.at("w").i.data(), fpoints.at("w").j.data(), fpoints.at("w").k.data(),                   // points to be forced
        fpoints.at("w").rot.data(),                                                                     // rotational matrix for local surface alignment
        fpoints.at("w").ip_u_i.data(), fpoints.at("w").ip_u_j.data(), fpoints.at("w").ip_u_k.data(), fpoints.at("u").c_idw_u.data(), // locations of the neighbouring u-points + weights
        fpoints.at("w").ip_v_i.data(), fpoints.at("w").ip_v_j.data(), fpoints.at("w").ip_v_k.data(), fpoints.at("u").c_idw_v.data(), // locations of the neighbouring v-points + weights
        fpoints.at("w").ip_w_i.data(), fpoints.at("w").ip_w_j.data(), fpoints.at("w").ip_w_k.data(), fpoints.at("u").c_idw_w.data(), // locations of the neighbouring w-points + weights
        fpoints.at("w").ip_s_i.data(), fpoints.at("w").ip_s_j.data(), fpoints.at("w").ip_s_k.data(), fpoints.at("u").c_idw_s.data(), // locations of the neighbouring s-points + weights
        fpoints.at("w").db.data(),                                                                      // distance nearest immersed boundary point to forcing point
        fpoints.at("w").di.data(),                                                                      // distance interpolation point to forcing point
        fpoints.at("w").z0b.data(),                                                                     // local roughness lengths of forcing points
        Boundary_type::Dirichlet_type,                                                                  // only allow no-slip conditions for momentum (for now..)
        fields.visc, fpoints.at("w").i.size(), n_idw_points,
        gd.icells, gd.ijcells,
        dt);

    set_ib_points(
        fields.mt.at("w")->fld.data(),
        fields.mp.at("w")->fld.data(),
        fpoints.at("w").ijk_ib.data(), TF(0.), 
        fpoints.at("w").ijk_ib.size(),
        dt);

    // not required here... CHECK, SvdL, 15-06-2023
    // boundary_cyclic.exec(fields.mt.at("u")->fld.data());
    // boundary_cyclic.exec(fields.mt.at("v")->fld.data());
    // boundary_cyclic.exec(fields.mt.at("w")->fld.data());

    for (auto &it : fields.sp) // SvdL, 21-05-2023: waarom staat hier een & voor it?
    {
        set_forcing_points_c(
            fields.st.at(it.first)->fld.data(),
            it.second->fld.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("w")->fld.data(),
            fpoints.at("s").sbot.at(it.first).data(),                                                       // value of boundary conditions to be enforced
            fpoints.at("s").i.data(), fpoints.at("s").j.data(), fpoints.at("s").k.data(),                   // points to be forced
            fpoints.at("s").rot.data(),                                                                     // rotational matrix for local surface alignment (although not used yet for scalars: DO calculate, needed for evisc)
            fpoints.at("s").ip_u_i.data(), fpoints.at("s").ip_u_j.data(), fpoints.at("s").ip_u_k.data(), fpoints.at("s").c_idw_u.data(), // locations of the neighbouring u-points + weights
            fpoints.at("s").ip_v_i.data(), fpoints.at("s").ip_v_j.data(), fpoints.at("s").ip_v_k.data(), fpoints.at("s").c_idw_v.data(), // locations of the neighbouring v-points + weights
            fpoints.at("s").ip_w_i.data(), fpoints.at("s").ip_w_j.data(), fpoints.at("s").ip_w_k.data(), fpoints.at("s").c_idw_w.data(), // locations of the neighbouring w-points + weights
            fpoints.at("s").ip_s_i.data(), fpoints.at("s").ip_s_j.data(), fpoints.at("s").ip_s_k.data(), fpoints.at("s").c_idw_s.data(), // locations of the neighbouring s-points + weights
            fpoints.at("s").db.data(),                                                                      // distance nearest immersed boundary point to forcing point
            fpoints.at("s").di.data(),                                                                      // distance interpolation point to forcing point
            fpoints.at("s").z0b.data(),                                                                     // local roughness lengths of forcing points (all scalars will have same for now..)
            sbc.at(it.first),                                                                                         // should contain Boundary_Type:: for all scalars (make variation between scalars possible?)
            fields.visc, fpoints.at("s").i.size(), n_idw_points,
            gd.icells, gd.ijcells,
            dt);

        // SvdL, 24-07-2024, backup function above: sbcbot,       // should contain Boundary_Type:: for all scalars (make variation between scalars possible?)
        // SvdL, 14-06-2023: the value to force inside still has to be set properly..
        set_ib_points(
            fields.st.at(it.first)->fld.data(),
            it.second->fld.data(),
            fpoints.at("s").ijk_ib.data(), TF(0.), 
            fpoints.at("s").ijk_ib.size(),
            dt);

        // not required here... CHECK, SvdL, 15-06-2023
        // boundary_cyclic.exec(it.second->fld.data());
        // boundary_cyclic.exec(fields.st.at(it.first)->fld.data());
    }

    // SvdL, 21-05-203: plans for much much later... allowing for IB conditions to be updated
}
#endif

// SvdL, 18-05-2023: fix later..
template <typename TF>
void Immersed_boundary<TF>::init(Input &inputin, Cross<TF> &cross)
{
    auto &gd = grid.get_grid_data();

    if (sw_ib == IB_type::Disabled)
        return;

    // // SvdL, 24-05-2023: for now set this to Dirichlet for all scalars.. eventually put this in the fpoints structure per scalar.
    // sbcbot = Boundary_type::Dirichlet_type;

    // Initialize structures for forcing points
    fpoints.emplace("u", Forcing_points<TF>());
    fpoints.emplace("v", Forcing_points<TF>());
    fpoints.emplace("w", Forcing_points<TF>());
    fpoints.emplace("s", Forcing_points<TF>()); //<< always initialize one for scalars (eddy diffusivity is at this location)

    // TENTATIVE SvdL, 23-07-2024: fix placement of different scalar fields (incl. temperature later). Currently, not expected to work for thermo
    // DEFINITELY NOT THE NICEST IMPLEMENTATION... I JUST ASSUME ALSO THAT ORDERING OF SCALAR FIELDS IS EQUAL TO READ IN OF SWBOTLIST
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

            ++n;
        }
    }
    
    // SvdL, 22-05-2023: niet geheel duidelijk wat dit nu doet..
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
    else if (sw_ib == IB_type::STL)
    {
        return;
    }
    else if (sw_ib == IB_type::User)
    {
        // Get grid and MPI information
        auto &gd = grid.get_grid_data();
        auto &mpi = master.get_MPI_data();

        const int ii = 1;
        const int jj = gd.icells;
        const int kk = gd.ijcells;

        // Offsets used in the the force points attribution
        const int mpi_offset_x = mpi.mpicoordx * gd.imax; // + gd.igc; //SvdL, 02-06-2023: pretty sure this should be a plus sign and without adding ghost cells, they are added below
        const int mpi_offset_y = mpi.mpicoordy * gd.jmax; // + gd.jgc;

        const int i0 = mpi_offset_x + gd.igc;             // SvdL, 02-06-2023: these represent the "beginning" and "end" values of i, j in the complete grid indexing (as if no MPI was present)
        const int i1 = mpi_offset_x + gd.igc + gd.imax;
        const int j0 = mpi_offset_y + gd.jgc;
        const int j1 = mpi_offset_y + gd.jgc + gd.jmax;

        const int rdim = 9;

        //         // For Debugging, SvdL 24 juli 2024
        // volatile int i = 0;
        // char hostname[256];
        // gethostname(hostname, sizeof(hostname));
        // printf("PID %d on %s ready for attach\n", getpid(), hostname);
        // fflush(stdout);
        // while (0 == i)
        //     sleep(5);

        // Get IBM settings from NetCDF input.
        Netcdf_group &init_group = input_nc.get_group("ibm");

        for (auto &it : fields.mp)
        {
            // 1. Make temporary vectors and read all data
            std::string n_fpoints_dim = "nfpoints_" + it.first;
            const int n_fpoints = init_group.get_dimension_size(n_fpoints_dim);

            std::vector<int> i_fpoints(n_fpoints);
            std::vector<int> j_fpoints(n_fpoints); 
            std::vector<int> k_fpoints(n_fpoints); 

            // std::vector<TF> xb_fpoints(n_fpoints);  
            // std::vector<TF> yb_fpoints(n_fpoints);  
            // std::vector<TF> zb_fpoints(n_fpoints);  

            std::vector<TF> xi_fpoints(n_fpoints);  
            std::vector<TF> yi_fpoints(n_fpoints);
            std::vector<TF> zi_fpoints(n_fpoints);

            std::vector<TF> db_fpoints(n_fpoints); 
            std::vector<TF> di_fpoints(n_fpoints); 
            std::vector<TF> z0b_fpoints(n_fpoints); 
            std::vector<TF> mbot_fpoints(n_fpoints);
            
            std::vector<TF> nor_fpoints(3 * n_fpoints);
           
            std::vector<int> start = {0};
            std::vector<int> count = {n_fpoints}; // SvdL, 01-06-2023: check if this total counter is correct

            init_group.get_variable(i_fpoints, "i_" + it.first, start, count);
            init_group.get_variable(j_fpoints, "j_" + it.first, start, count);
            init_group.get_variable(k_fpoints, "k_" + it.first, start, count);

            // init_group.get_variable(xb_fpoints, "xb_"+ it.first, start, count);
            // init_group.get_variable(yb_fpoints, "yb_"+ it.first, start, count);
            // init_group.get_variable(zb_fpoints, "zb_"+ it.first, start, count);

            init_group.get_variable(xi_fpoints, "xi_"+ it.first, start, count);
            init_group.get_variable(yi_fpoints, "yi_"+ it.first, start, count);
            init_group.get_variable(zi_fpoints, "zi_"+ it.first, start, count);

            init_group.get_variable(db_fpoints, "db_" + it.first, start, count);
            init_group.get_variable(di_fpoints, "di_" + it.first, start, count);

            init_group.get_variable(z0b_fpoints, "z0b_" + it.first, start, count);
            init_group.get_variable(mbot_fpoints, "mbot_" + it.first, start, count);

            count = {3 * n_fpoints};   // SvdL, 01-06-2023: check if this total counter is correct, and if needed..
            init_group.get_variable(nor_fpoints, "nor_" + it.first, start, count);

            // count = {rdim * n_fpoints}; // SvdL, 01-06-2023: check if this total counter is correct
            // init_group.get_variable(rot_fpoints, "rot_" + it.first, start, count);

            // 2. For all entries, check whether on current MPI task
            for (int nn=0; nn<n_fpoints; ++nn)
            {
                int i = i_fpoints[nn] + gd.igc; // SvdL, 02-06-2023: these only have local scope
                int j = j_fpoints[nn] + gd.jgc; // do add ghost cells, as also for serial job all i,j,k indices will have been shifted by nr of ghost cells.
                int k = k_fpoints[nn] + gd.kgc; // could remove ghost cell here AND above at i0, i1, j0 and j1 declerations. That would cancel out.. << NO PROBABLY WRONG TO DO SO (SvdL, 25-07-2024)

                // Check if grid point on this MPI task
                if (i >= i0 && i < i1 && j >= j0 && j < j1)
                {
                    i -= mpi_offset_x; // if check is passed, get i,j indices as would be at current MPI task
                    j -= mpi_offset_y;

                    const int ijk = i + j*jj + k*kk;

                    // 3. Start pushing back IBM data to correct Forcing_points structure 
                    fpoints.at(it.first).ijk_fp.push_back(ijk);

                    fpoints.at(it.first).i.push_back(i);
                    fpoints.at(it.first).j.push_back(j);
                    fpoints.at(it.first).k.push_back(k);

                    // fpoints.at(it.first).xb.push_back(xb_fpoints[nn]);
                    // fpoints.at(it.first).yb.push_back(yb_fpoints[nn]);
                    // fpoints.at(it.first).zb.push_back(zb_fpoints[nn]);

                    fpoints.at(it.first).xi.push_back(xi_fpoints[nn]);
                    fpoints.at(it.first).yi.push_back(yi_fpoints[nn]);
                    fpoints.at(it.first).zi.push_back(zi_fpoints[nn]);

                    fpoints.at(it.first).db.push_back(db_fpoints[nn]);
                    fpoints.at(it.first).di.push_back(di_fpoints[nn]);

                    fpoints.at(it.first).z0b.push_back(z0b_fpoints[nn]);
                    fpoints.at(it.first).mbot.push_back(mbot_fpoints[nn]);

                    fpoints.at(it.first).nor.push_back(nor_fpoints[3*nn]);   // SvdL, 02-06-2023: check these..
                    fpoints.at(it.first).nor.push_back(nor_fpoints[3*nn+1]);
                    fpoints.at(it.first).nor.push_back(nor_fpoints[3*nn+2]);

                    // also still include components of rotation vector..? Or internally calculate in microhh..?
                    // SvdL, 02-06-2023: stupid note here, DO replace nfoints values everywhere with some get_size operation..
                }
            }
            
            // 4. now that amount of forcing points at MPI task are known, preset containers for interpolation points and weights, and rotation matrix
            const int fpoints_on_mpi_task = fpoints.at(it.first).i.size();

            fpoints.at(it.first).ip_u_i.resize(fpoints_on_mpi_task*n_idw_points);
            fpoints.at(it.first).ip_u_j.resize(fpoints_on_mpi_task*n_idw_points);
            fpoints.at(it.first).ip_u_k.resize(fpoints_on_mpi_task*n_idw_points);
            fpoints.at(it.first).ip_v_i.resize(fpoints_on_mpi_task*n_idw_points);
            fpoints.at(it.first).ip_v_j.resize(fpoints_on_mpi_task*n_idw_points);
            fpoints.at(it.first).ip_v_k.resize(fpoints_on_mpi_task*n_idw_points);
            fpoints.at(it.first).ip_w_i.resize(fpoints_on_mpi_task*n_idw_points);
            fpoints.at(it.first).ip_w_j.resize(fpoints_on_mpi_task*n_idw_points);
            fpoints.at(it.first).ip_w_k.resize(fpoints_on_mpi_task*n_idw_points);
            fpoints.at(it.first).ip_s_i.resize(fpoints_on_mpi_task*n_idw_points);
            fpoints.at(it.first).ip_s_j.resize(fpoints_on_mpi_task*n_idw_points);
            fpoints.at(it.first).ip_s_k.resize(fpoints_on_mpi_task*n_idw_points);

            fpoints.at(it.first).c_idw_u.resize(fpoints_on_mpi_task*n_idw_points);
            fpoints.at(it.first).c_idw_v.resize(fpoints_on_mpi_task*n_idw_points);
            fpoints.at(it.first).c_idw_w.resize(fpoints_on_mpi_task*n_idw_points);
            fpoints.at(it.first).c_idw_s.resize(fpoints_on_mpi_task*n_idw_points);

            fpoints.at(it.first).rot.resize(fpoints_on_mpi_task*rdim);
        }

        // Split loading in of scalar points in two parts - first part is always needed for eddy viscosity..
        // (can be later changed for perfect DNS mode..)
        // 1. Make temporary vectors and read all data
        std::string n_fpoints_dim = "nfpoints_s";
        const int n_fpoints = init_group.get_dimension_size(n_fpoints_dim);

        std::vector<int> i_fpoints(n_fpoints);
        std::vector<int> j_fpoints(n_fpoints); 
        std::vector<int> k_fpoints(n_fpoints); 

        // std::vector<TF> xb_fpoints(n_fpoints);  
        // std::vector<TF> yb_fpoints(n_fpoints);  
        // std::vector<TF> zb_fpoints(n_fpoints);  

        std::vector<TF> xi_fpoints(n_fpoints);  
        std::vector<TF> yi_fpoints(n_fpoints);
        std::vector<TF> zi_fpoints(n_fpoints);

        std::vector<TF> db_fpoints(n_fpoints); 
        std::vector<TF> di_fpoints(n_fpoints); 
        std::vector<TF> z0b_fpoints(n_fpoints); 
        std::vector<TF> mbot_fpoints(n_fpoints);
        
        std::vector<TF> nor_fpoints(3 * n_fpoints);
        
        std::vector<int> start = {0};
        std::vector<int> count = {n_fpoints}; // SvdL, 01-06-2023: check if this total counter is correct

        init_group.get_variable(i_fpoints, "i_s", start, count);
        init_group.get_variable(j_fpoints, "j_s", start, count);
        init_group.get_variable(k_fpoints, "k_s", start, count);

        // init_group.get_variable(xb_fpoints, "xb_s", start, count);
        // init_group.get_variable(yb_fpoints, "yb_s", start, count);
        // init_group.get_variable(zb_fpoints, "zb_s", start, count);

        init_group.get_variable(xi_fpoints, "xi_s", start, count);
        init_group.get_variable(yi_fpoints, "yi_s", start, count);
        init_group.get_variable(zi_fpoints, "zi_s", start, count);

        init_group.get_variable(db_fpoints, "db_s", start, count);
        init_group.get_variable(di_fpoints, "di_s", start, count);

        init_group.get_variable(z0b_fpoints, "z0b_s", start, count);
        init_group.get_variable(mbot_fpoints, "mbot_s", start, count);

        count = {3 * n_fpoints};   // SvdL, 01-06-2023: check if this total counter is correct, and if needed..
        init_group.get_variable(nor_fpoints, "nor_s", start, count);

        // 2. For all entries, check whether on current MPI task
        for (int nn=0; nn<n_fpoints; ++nn)
        {
            int i = i_fpoints[nn] + gd.igc; // SvdL, 02-06-2023: these only have local scope
            int j = j_fpoints[nn] + gd.jgc; // do add ghost cells, as also for serial job all i,j,k indices will have been shifted by nr of ghost cells.
            int k = k_fpoints[nn] + gd.kgc; // could remove ghost cell here AND above at i0, i1, j0 and j1 declerations. That would cancel out..

            // Check if grid point on this MPI task
            if (i >= i0 && i < i1 && j >= j0 && j < j1)
            {
                i -= mpi_offset_x; // if check is passed, get i,j indices as would be at current MPI task
                j -= mpi_offset_y;

                const int ijk = i + j*jj + k*kk;

                // 3. Start pushing back IBM data to correct Forcing_points structure 
                fpoints.at("s").ijk_fp.push_back(ijk);

                fpoints.at("s").i.push_back(i);
                fpoints.at("s").j.push_back(j);
                fpoints.at("s").k.push_back(k);

                // fpoints.at("s").xb.push_back(xb_fpoints[nn]);
                // fpoints.at("s").yb.push_back(yb_fpoints[nn]);
                // fpoints.at("s").zb.push_back(zb_fpoints[nn]);

                fpoints.at("s").xi.push_back(xi_fpoints[nn]);
                fpoints.at("s").yi.push_back(yi_fpoints[nn]);
                fpoints.at("s").zi.push_back(zi_fpoints[nn]);

                fpoints.at("s").db.push_back(db_fpoints[nn]);
                fpoints.at("s").di.push_back(di_fpoints[nn]);

                fpoints.at("s").z0b.push_back(z0b_fpoints[nn]);
                // fpoints.at("s").mbot.push_back(mbot_fpoints[nn]); // SvdL, 05-06-2023: not relevant for scalars

                fpoints.at("s").nor.push_back(nor_fpoints[3*nn]);    
                fpoints.at("s").nor.push_back(nor_fpoints[3*nn+1]);
                fpoints.at("s").nor.push_back(nor_fpoints[3*nn+2]);

            }
        }

        // 4. now that amount of forcing points at MPI task are known, preset containers for interpolation points and weights, and rotation matrix
        const int fpoints_on_mpi_task = fpoints.at("s").i.size();

        fpoints.at("s").ip_u_i.resize(fpoints_on_mpi_task*n_idw_points);
        fpoints.at("s").ip_u_j.resize(fpoints_on_mpi_task*n_idw_points);
        fpoints.at("s").ip_u_k.resize(fpoints_on_mpi_task*n_idw_points);
        fpoints.at("s").ip_v_i.resize(fpoints_on_mpi_task*n_idw_points);
        fpoints.at("s").ip_v_j.resize(fpoints_on_mpi_task*n_idw_points);
        fpoints.at("s").ip_v_k.resize(fpoints_on_mpi_task*n_idw_points);
        fpoints.at("s").ip_w_i.resize(fpoints_on_mpi_task*n_idw_points);
        fpoints.at("s").ip_w_j.resize(fpoints_on_mpi_task*n_idw_points);
        fpoints.at("s").ip_w_k.resize(fpoints_on_mpi_task*n_idw_points);
        fpoints.at("s").ip_s_i.resize(fpoints_on_mpi_task*n_idw_points);
        fpoints.at("s").ip_s_j.resize(fpoints_on_mpi_task*n_idw_points);
        fpoints.at("s").ip_s_k.resize(fpoints_on_mpi_task*n_idw_points);

        fpoints.at("s").c_idw_u.resize(fpoints_on_mpi_task*n_idw_points);
        fpoints.at("s").c_idw_v.resize(fpoints_on_mpi_task*n_idw_points);
        fpoints.at("s").c_idw_w.resize(fpoints_on_mpi_task*n_idw_points);
        fpoints.at("s").c_idw_s.resize(fpoints_on_mpi_task*n_idw_points);

        fpoints.at("s").rot.resize(fpoints_on_mpi_task*rdim);

        // this second part is only required if there are actually scalars (including temperature..)
        if (fields.sp.size() > 0)
        {
            // Set scalar boundary conditions per scalar
            // .. unfortunately, this requires looping over nn again. Can this be improved?
            std::vector<TF> sbot_fpoints(n_fpoints);
            start = {0};
            count = {n_fpoints}; // SvdL, 23-05-2023: check if this total counter is correct
        
            for (auto &scalar : fields.sp)
            {
                init_group.get_variable(sbot_fpoints, "sbot_" + scalar.first, start, count); // SvdL, 24-05-2023: ik verwacht dat deze assignment niet mag..

                for (int nn=0; nn<n_fpoints; ++nn)
                {
                    int i = i_fpoints[nn] + gd.igc; // SvdL, 02-06-2023: these only have local scope
                    int j = j_fpoints[nn] + gd.jgc; // do add ghost cells, as also for serial job all i,j,k indices will have been shifted by nr of ghost cells.
                    int k = k_fpoints[nn] + gd.kgc; // could remove ghost cell here AND above at i0, i1, j0 and j1 declerations. That would cancel out..

                    // Check if grid point on this MPI task
                    if (i >= i0 && i < i1 && j >= j0 && j < j1)
                    {
                        i -= mpi_offset_x; // if check is passed, get i,j indices at would be at current MPI task
                        j -= mpi_offset_y;

                        const int ijk = i + j*jj + k*kk;

                        fpoints.at("s").sbot.at(scalar.first).push_back(sbot_fpoints[nn]);
                        
                    }
                }
            }
        }

        // Svd:, 14-06-2023: separately read in the indices of the ib points
        // is probably easier to just read a full mask (requires however additional memory.. 4 fields)
        for (auto &it : fields.mp)
        { 
            std::string n_ib_dim = "ibpoints_" + it.first;
            const int n_ib = init_group.get_dimension_size(n_ib_dim);
        
            std::vector<int> i_ibpoints(n_ib); 
            std::vector<int> j_ibpoints(n_ib); 
            std::vector<int> k_ibpoints(n_ib); 

            std::vector<int> start = {0};
            std::vector<int> count = {n_ib}; // SvdL, 01-06-2023: check if this total counter is correct

            init_group.get_variable(i_ibpoints, "i_ib_"+ it.first, start, count);
            init_group.get_variable(j_ibpoints, "j_ib_"+ it.first, start, count);
            init_group.get_variable(k_ibpoints, "k_ib_"+ it.first, start, count);

            // 2. For all entries, check whether on current MPI task
            for (int nn=0; nn<n_ib; ++nn)
            {
                int i = i_ibpoints[nn] + gd.igc; // SvdL, 02-06-2023: these only have local scope
                int j = j_ibpoints[nn] + gd.jgc; // do add ghost cells, as also for serial job all i,j,k indices will have been shifted by nr of ghost cells.
                int k = k_ibpoints[nn] + gd.kgc; // could remove ghost cell here AND above at i0, i1, j0 and j1 declerations. That would cancel out..

                // Check if grid point on this MPI task
                if (i >= i0 && i < i1 && j >= j0 && j < j1)
                {
                    i -= mpi_offset_x; // if check is passed, get i,j indices as would be at current MPI task
                    j -= mpi_offset_y;

                    const int ijk = i + j*jj + k*kk;
                    fpoints.at(it.first).ijk_ib.push_back(ijk);
            
                }
            }        
        }
        
        // Finally, repeat this again for the scalar position
        std::string n_ib_dim = "ibpoints_s";
        const int n_ib = init_group.get_dimension_size(n_ib_dim);

        std::vector<int> i_ibpoints(n_ib); 
        std::vector<int> j_ibpoints(n_ib); 
        std::vector<int> k_ibpoints(n_ib); 

        start = {0};
        count = {n_ib}; // SvdL, 01-06-2023: check if this total counter is correct

        init_group.get_variable(i_ibpoints, "i_ib_s", start, count);
        init_group.get_variable(j_ibpoints, "j_ib_s", start, count);
        init_group.get_variable(k_ibpoints, "k_ib_s", start, count);

        // 2. For all entries, check whether on current MPI task
        for (int nn=0; nn<n_ib; ++nn)
        {
            int i = i_ibpoints[nn] + gd.igc; // SvdL, 02-06-2023: these only have local scope
            int j = j_ibpoints[nn] + gd.jgc; // do add ghost cells, as also for serial job all i,j,k indices will have been shifted by nr of ghost cells.
            int k = k_ibpoints[nn] + gd.kgc; // could remove ghost cell here AND above at i0, i1, j0 and j1 declerations. That would cancel out..

            // Check if grid point on this MPI task
            if (i >= i0 && i < i1 && j >= j0 && j < j1)
            {
                i -= mpi_offset_x; // if check is passed, get i,j indices as would be at current MPI task
                j -= mpi_offset_y;

                const int ijk = i + j*jj + k*kk;

                fpoints.at("s").ijk_ib.push_back(ijk);

            }
        }        

        // Init the toolbox classes.
        boundary_cyclic.init();
        
        for (auto &it : fields.mp)
        {
            setup_interpolation(
                fpoints.at(it.first).xi.data(), fpoints.at(it.first).yi.data(), fpoints.at(it.first).zi.data(),
                fpoints.at(it.first).ip_u_i.data(), fpoints.at(it.first).ip_u_j.data(), fpoints.at(it.first).ip_u_k.data(), fpoints.at(it.first).c_idw_u.data(), 
                fpoints.at(it.first).ip_v_i.data(), fpoints.at(it.first).ip_v_j.data(), fpoints.at(it.first).ip_v_k.data(), fpoints.at(it.first).c_idw_v.data(), 
                fpoints.at(it.first).ip_w_i.data(), fpoints.at(it.first).ip_w_j.data(), fpoints.at(it.first).ip_w_k.data(), fpoints.at(it.first).c_idw_w.data(), 
                fpoints.at(it.first).ip_s_i.data(), fpoints.at(it.first).ip_s_j.data(), fpoints.at(it.first).ip_s_k.data(), fpoints.at(it.first).c_idw_s.data(), 
                fpoints.at("u").ijk_fp, fpoints.at("u").ijk_ib,  //<< this passing just has to do with the structure in which the data is saved
                fpoints.at("v").ijk_fp, fpoints.at("v").ijk_ib,
                fpoints.at("w").ijk_fp, fpoints.at("w").ijk_ib,
                fpoints.at("s").ijk_fp, fpoints.at("s").ijk_ib,
                fpoints.at(it.first).i.size(), n_idw_points_min, n_idw_points, 
                gd.x, gd.y, gd.z, gd.xh, gd.yh, gd.zh, gd.dx, gd.dy, gd.dz,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend, gd.jend, gd.kend,
                gd.icells, gd.ijcells, 
                mpi_offset_x, mpi_offset_y);

            calculate_rotation_matrix(
                fpoints.at(it.first).rot.data(), 
                fpoints.at(it.first).nor.data(), 
                fpoints.at(it.first).i.size());

        }

        setup_interpolation(
            fpoints.at("s").xi.data(), fpoints.at("s").yi.data(), fpoints.at("s").zi.data(),
            fpoints.at("s").ip_u_i.data(), fpoints.at("s").ip_u_j.data(), fpoints.at("s").ip_u_k.data(), fpoints.at("s").c_idw_u.data(), 
            fpoints.at("s").ip_v_i.data(), fpoints.at("s").ip_v_j.data(), fpoints.at("s").ip_v_k.data(), fpoints.at("s").c_idw_v.data(), 
            fpoints.at("s").ip_w_i.data(), fpoints.at("s").ip_w_j.data(), fpoints.at("s").ip_w_k.data(), fpoints.at("s").c_idw_w.data(), 
            fpoints.at("s").ip_s_i.data(), fpoints.at("s").ip_s_j.data(), fpoints.at("s").ip_s_k.data(), fpoints.at("s").c_idw_s.data(), 
            fpoints.at("u").ijk_fp, fpoints.at("u").ijk_ib,  //<< this passing just has to do with the structure in which the data is saved
            fpoints.at("v").ijk_fp, fpoints.at("v").ijk_ib,
            fpoints.at("w").ijk_fp, fpoints.at("w").ijk_ib,
            fpoints.at("s").ijk_fp, fpoints.at("s").ijk_ib,
            fpoints.at("s").i.size(), n_idw_points_min, n_idw_points, 
            gd.x, gd.y, gd.z, gd.xh, gd.yh, gd.zh, gd.dx, gd.dy, gd.dz,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend, gd.jend, gd.kend,
            gd.icells, gd.ijcells,
            mpi_offset_x, mpi_offset_y);

        calculate_rotation_matrix(
            fpoints.at("s").rot.data(), 
            fpoints.at("s").nor.data(), 
            fpoints.at("s").i.size());

    }

    // SvdL, 21-05-2023: output some important notes.
    master.print_message("SvdL: Realize that statistics between 1st and 2nd layer are likely meaningless at this stage! (Because IB-values are only set after routines statistics are calculated\n");
    master.print_message("SvdL: further, the statistics at the forcing points need to be altered to the IBM forcing, with all other tendencies set to zero here.\n");
    //     // Print some statistics (number of ghost cells)
    //     print_statistics(ghost.at("u").i, std::string("u"), master);
}

// // SvdL, 18-05-2023: fix later.. for now make it BS-function
template <typename TF>
bool Immersed_boundary<TF>::has_mask(std::string mask_name)
{
    // if (std::find(available_masks.begin(), available_masks.end(), mask_name) != available_masks.end())
    //     return true;
    // else
    //     return false;
    return false;
}

// // SvdL, 18-05-2023: fix later.. or now make it BS-function
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

// SvdL, 18-05-2023: fix later..
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

/* IMPLEMENTATION NOTES - STATE 15-06-2023

This changed IBM implementation is based on Roman et al. (200x), in which the first layer of fluid points outside the IB are locally forced to comply with 
the logarithmic law of the wall. 

The ib forcings are calculated using the so-called auxiliary velocity and scalars, meaning that we first calculate what would be the fluid velocities without ibm 
and without pressure (fluctuation) gradient forcing. These intermediate velocities are then used to determine the velocities and scalars at the forcing points. 
These new values and the old values at the forcing points are finally used to overwrite the tendencies at the forcing points itself.

Taking this approach, systematic errors in the long run are hopefully suppressed. Using the velocity at the current timestep itself, the velocities at the forcing
points won't "match" with those at the fluid points at the next timestep. We use the knowlegde that the auxiliary velocity will be a "better" predictor of the
velocities in the fluid at the next timestep than the current values.

We DO have to apply the ibm tendencies before the pressure solver, as we DO want to still globally enforce a divergence free flow.

OLD NOTE, controleer nog..
NOTE: DIT IS EIGENLIJK ILL-DEFINED... JUIST OMDAT DE INDEX 1 KAN VERSCHILLEN (MAX) MOET JE DE SEARCH REGION UITBREIDEN OM HIER REKENING MEE TE HOUDEN.. 
DIT KAN ER VOOR ZORGEN DAT JE DAARDOOR JUIST PUNTEN TE VER WEG MEE GAAT NEMEN, OF NIET???
VERDERE NOTE: METHODE NEEMT IMPLICIET OOK AAN DAT JE TUSSEN 2DE EN 3DE GRIDPUNT VAN DE IB GENOEG RESOLUTIE, DAN WEL AFNAME GRADIENT HEBT Z.D.D. LOGGRADIENT RESOLVED WORDT!!

*/

template class Immersed_boundary<double>;
template class Immersed_boundary<float>;
