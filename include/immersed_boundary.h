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

enum class IB_type {Disabled, SDF};

// Forcing points info on staggered grid (for SDF implementation, IB_type=SDF)
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

    // Combined grid indices ijk for forcing point locations
    std::vector<int> ijk;

    // Nearest location of IB to forcing cell:
    std::vector<TF> xb;     // size = number of forcing points
    std::vector<TF> yb;
    std::vector<TF> zb;

    // Location of interpolation point on normal from IB to forcing cell:
    std::vector<TF> xi;     // size = number of forcing points
    std::vector<TF> yi;
    std::vector<TF> zi;

    // Local normal vectors and surface rotation matrix:
    // std::vector<TF> nor;   // size = number of forcing points x 3
    std::vector<TF> rot;   // size = number of forcing points x 9

    // Distance forcing points to immersed boundary (dist_b) and interpolation point (dist_i)
    std::vector<TF> dist_b; 
    std::vector<TF> dist_i; 

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
    
    // //
    // // GPU 
    // //

    // Indices of the points to be forced:
    int* i_g;     // size = number of forcing points
    int* j_g;
    int* k_g;

    // // Combined grid indices ijk for forcing point locations
    // int* ijk_g;

    // // Nearest location of IB to forcing cell:
    // TF* xb_g;     // size = number of forcing points
    // TF* yb_g;
    // TF* zb_g;

    // // Location of interpolation point on normal from IB to forcing cell:
    // TF* xi_g;     // size = number of forcing points
    // TF* yi_g;
    // TF* zi_g;

    // Local normal vectors and surface rotation matrix:
    // std::vector<TF> nor;   // size = number of forcing points x 3
    TF* rot_g;   // size = number of forcing points x 9

    // Distance forcing points to immersed boundary (dist_b) and interpolation point (dist_i)
    TF* dist_b_g; 
    TF* dist_i_g; 

    // Aerodynamic roughness at immersed boundary / forcing point (can be different for momentum and scalars)
    TF* z0b_g;

    // Points outside IB used for IDW interpolation: //<  SvdL, 18-05-2023: remove probably the inversed distance weighting??
    int* ip_u_i_g;  // size = number of forcing points x n_idw_points (SvdL, check this comment)
    int* ip_u_j_g;
    int* ip_u_k_g;
    int* ip_v_i_g;  // size = number of forcing points x n_idw_points (SvdL, check this comment)
    int* ip_v_j_g;
    int* ip_v_k_g;
    int* ip_w_i_g;  // size = number of forcing points x n_idw_points (SvdL, check this comment)
    int* ip_w_j_g;
    int* ip_w_k_g;
    int* ip_s_i_g;  // size = number of forcing points x n_idw_points (SvdL, check this comment)
    int* ip_s_j_g;
    int* ip_s_k_g;

    // Interpolation coefficients of neighbours to interpolation point: 4 vectors due to staggered grid nature.
    TF* c_idw_u_g;     // size = number of ghost points x n_idw_points
    TF* c_idw_v_g;     // size = number of ghost points x n_idw_points
    TF* c_idw_w_g;     // size = number of ghost points x n_idw_points
    TF* c_idw_s_g;     // size = number of ghost points x n_idw_points

    // Spatially varying scalar (and momentum..) boundary conditions
    std::map<std::string, TF*> sbot_g;
    TF* mbot_g;
    
};

template<typename TF>
struct IB_points
{
    int n_ibpoints;             // number of ib points (i.e., within object)

    //
    // CPU
    //

    // Indices of the ib points:
    std::vector<int> i;        // size = number of ib points
    std::vector<int> j;
    std::vector<int> k;
    std::vector<int> ijk;

    // //
    // // GPU 
    // //

    // Indices of the points to be forced:
    int* i_g;     // size = number of ib points
    int* j_g;
    int* k_g;
    int* ijk_g;

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
        void exec(double);

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

        // All forcing points and ib properties
        std::map<std::string, Forcing_points<TF>> fpoints;
        std::map<std::string, IB_points<TF>> ibpoints;

        // Statistics
        std::vector<std::string> available_masks;

        // Cross-sections
        std::vector<std::string> crosslist;

        // Additional parameters for ib
        TF z0bound;
};

#endif
