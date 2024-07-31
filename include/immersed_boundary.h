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

#ifndef IMMERSED_BOUNDARY_H
#define IMMERSED_BOUNDARY_H

#include "boundary.h"
#include "boundary_cyclic.h"

class Master;
class Input;
class Netcdf_handle;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Cross;
// template<typename> class Diff;

enum class IB_type {Disabled, STL, User};

// Forcing points info on staggered grid
template<typename TF>
struct Forcing_points
{
    int n_fpoints;

    //
    // CPU
    //

    // Indices of the points to be forced:
    std::vector<int> i;     // size = number of forcing points
    std::vector<int> j;
    std::vector<int> k;

    // Nearest location of IB to forcing cell:
    std::vector<TF> xb;     // size = number of forcing points
    std::vector<TF> yb;
    std::vector<TF> zb;

    // Location of interpolation point on normal from IB to forcing cell:
    std::vector<TF> xi;     // size = number of forcing points
    std::vector<TF> yi;
    std::vector<TF> zi;

    // Local normal vectors and surface rotation matrix:
    std::vector<TF> nor;   // size = number of forcing points x 3
    std::vector<TF> rot;   // size = number of forcing points x 9

    // Distance forcing points to immersed boundary (db) and interpolation point (di)
    std::vector<TF> db; 
    std::vector<TF> di; 

    // Aerodynamic roughness at immersed boundary / forcing point (can be different for momentum and scalars)
    std::vector<TF> z0b;

    // Points outside IB used for IDW interpolation: //<  SvdL, 18-05-2023: remove probably the inversed distance weighting??
    std::vector<int> ip_u_i;  // size = number of forcing points x n_idw_points (SvdL, check this comment)
    std::vector<int> ip_u_j;
    std::vector<int> ip_u_k;
    std::vector<int> ip_v_i;  // size = number of forcing points x n_idw_points (SvdL, check this comment)
    std::vector<int> ip_v_j;
    std::vector<int> ip_v_k;
    std::vector<int> ip_w_i;  // size = number of forcing points x n_idw_points (SvdL, check this comment)
    std::vector<int> ip_w_j;
    std::vector<int> ip_w_k;
    std::vector<int> ip_s_i;  // size = number of forcing points x n_idw_points (SvdL, check this comment)
    std::vector<int> ip_s_j;
    std::vector<int> ip_s_k;

    // Interpolation coefficients of neighbours to interpolation point: 4 vectors due to staggered grid nature.
    std::vector<TF> c_idw_u;     // size = number of ghost points x n_idw_points
    std::vector<TF> c_idw_v;     // size = number of ghost points x n_idw_points
    std::vector<TF> c_idw_w;     // size = number of ghost points x n_idw_points
    std::vector<TF> c_idw_s;     // size = number of ghost points x n_idw_points

    // Spatially varying scalar (and momentum..) boundary conditions
    std::map<std::string, std::vector<TF>> sbot;
    std::vector<TF> mbot;
    
    // Combined grid indices ijk for forcing point locations and locations inside the ib
    std::vector<int> ijk_fp;
    std::vector<int> ijk_ib;

    // // SvdL, 21-05-2023: for now restrict to getting it working for CPU...
    // //
    // // GPU 
    // //

    // // Indices of IB ghost points:
    // int* i_g;
    // int* j_g;
    // int* k_g;

    // // Nearest location of IB to ghost cell:
    // TF* xb_g;
    // TF* yb_g;
    // TF* zb_g;

    // // Location of interpolation point outside IB:
    // TF* xi_g;
    // TF* yi_g;
    // TF* zi_g;

    // TF* di_g;  // Distance ghost cell to interpolation point

    // // Points outside IB used for IDW interpolation:
    // int* ip_i_g;
    // int* ip_j_g;
    // int* ip_k_g;
    // TF* ip_d_g;

    // // Interpolation coefficients
    // TF* c_idw_g;
    // TF* c_idw_sum_g;

    // // Spatially varying scalar (and momentum..) boundary conditions
    // std::map<std::string, TF*> sbot_g;
    // TF* mbot_g;
};

// Convenience struct to simplify sorting
template<typename TF>
struct Neighbour
{
    int i;
    int j;
    int k;
    TF distance;
};

template<typename TF>
class Immersed_boundary
{
    public:
        Immersed_boundary(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Immersed_boundary();

        void init(Input&, Cross<TF>&);
        void create(Input&, Netcdf_handle&);

        void exec_viscosity();
        void exec(const double);

        void exec_cross(Cross<TF>&, unsigned long);

        bool has_mask(std::string);
        void get_mask(Stats<TF>&, std::string);

        void prepare_device();
        void clear_device();

        IB_type get_switch() const { return sw_ib; }


    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        // Diff<TF>& diff;
        Field3d_io<TF> field3d_io;
        Boundary_cyclic<TF> boundary_cyclic;

        IB_type sw_ib;

        int n_idw_points;       // Number of interpolation points in IDW interpolation
        int n_idw_points_min;   // Minimum number of interpolation points in IDW interpolation

        // Boundary type for scalars
        // SvdL, 20240724: for now, keep constant for all buildings points (later to be included in map sbot above)
        Boundary_type sbcbot;
        std::vector<std::string> swbotlist; //<< list containing boundary type per scalar (later incl. temperature separately!)
        std::map<std::string, Boundary_type> sbc; //<< not implemented yet
        
        // std::vector<std::string> sbot_spatial_list; //<< not implemented yet

        // All forcing points properties
        std::map<std::string, Forcing_points<TF>> fpoints;

        // Statistics
        std::vector<std::string> available_masks;

        // Cross-sections
        std::vector<std::string> crosslist;

        // SvdL, 23-05-2023: I assume this can be a private member.. anyway still implement
        // void process_ibm_input(Input&, Netcdf_handle&); 
};

#endif
