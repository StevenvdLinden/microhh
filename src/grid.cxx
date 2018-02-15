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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "master.h"
#include "grid.h"
#include "input.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "model.h"

#include <iostream> // REMOVE ME BvS

/**
 * This function constructs the grid class.
 * @param modelin Pointer to the model class.
 * @param inputin Pointer to the input class.
 */
Grid::Grid(Model *modelin, Input *inputin)
{
    master = modelin->master;

    mpitypes  = false;
    fftwplan  = false;

    // Initialize the pointers to zero.
    x  = 0;
    xh = 0;
    y  = 0;
    yh = 0;
    z  = 0;
    zh = 0;

    dz    = 0;
    dzh   = 0;
    dzi   = 0;
    dzhi  = 0;
    dzi4  = 0;
    dzhi4 = 0;

    z_g     = 0;
    zh_g    = 0;
    dz_g    = 0;
    dzh_g   = 0;
    dzi_g   = 0;
    dzhi_g  = 0;
    dzi4_g  = 0;
    dzhi4_g = 0;

    fftini  = 0;
    fftouti = 0;
    fftinj  = 0;
    fftoutj = 0;

    int nerror = 0;
    nerror += inputin->get_item(&xsize, "grid", "xsize", "");
    nerror += inputin->get_item(&ysize, "grid", "ysize", "");
    nerror += inputin->get_item(&zsize, "grid", "zsize", "");

    nerror += inputin->get_item(&itot, "grid", "itot", "");
    nerror += inputin->get_item(&jtot, "grid", "jtot", "");
    nerror += inputin->get_item(&ktot, "grid", "ktot", "");

    nerror += inputin->get_item(&utrans, "grid", "utrans", "", 0.);
    nerror += inputin->get_item(&vtrans, "grid", "vtrans", "", 0.);

    nerror += inputin->get_item(&swspatialorder, "grid", "swspatialorder", "");

    if (nerror)
        throw 1;

    if (!(swspatialorder == "2" || swspatialorder == "4"))
    {
        master->print_error("\"%s\" is an illegal value for swspatialorder\n", swspatialorder.c_str());
        throw 1;
    }
    // 2nd order scheme requires only 1 ghost cell
    if (swspatialorder == "2")
    {
        igc = 1;
        jgc = 1;
        kgc = 1;
    }
    // 4th order scheme requires 3 ghost cells
    else if (swspatialorder == "4")
    {
        igc = 3;
        jgc = 3;
        kgc = 3;
    }
}

/**
 * This function destructs the grid class.
 */
Grid::~Grid()
{
    if (fftwplan)
    {
        fftw_destroy_plan(iplanf);
        fftw_destroy_plan(iplanb);
        fftw_destroy_plan(jplanf);
        fftw_destroy_plan(jplanb);
    }

    delete[] x;
    delete[] xh;
    delete[] y;
    delete[] yh;
    delete[] z;
    delete[] zh;
    delete[] dz;
    delete[] dzh;
    delete[] dzi;
    delete[] dzhi;
    delete[] dzi4;
    delete[] dzhi4;

    fftw_free(fftini);
    fftw_free(fftouti);
    fftw_free(fftinj);
    fftw_free(fftoutj);

    fftw_cleanup();

#ifdef USECUDA
    clear_device();
#endif

    exit_mpi();
}

/**
 * This function allocates the dynamic arrays in the field class
 * variables and calculates the derived grid indices and dimensions.
 */
void Grid::init()
{
    // Check whether the grid fits the processor configuration.
    if (itot % master->npx != 0)
    {
        master->print_error("itot = %d is not a multiple of npx = %d\n", itot, master->npx);
        throw 1;
    }
    if (itot % master->npy != 0)
    {
        master->print_error("itot = %d is not a multiple of npy = %d\n", itot, master->npy);
        throw 1;
    }
    // Check this one only when npy > 1, since the transpose in that direction only happens then.
    if (jtot % master->npx != 0)
    {
        master->print_error("jtot = %d is not a multiple of npx = %d\n", jtot, master->npx);
        throw 1;
    }
    if (jtot % master->npy != 0)
    {
        master->print_error("jtot = %d is not a multiple of npy = %d\n", jtot, master->npy);
        throw 1;
    }
    if (ktot % master->npx != 0)
    {
        master->print_error("ERROR ktot = %d is not a multiple of npx = %d\n", ktot, master->npx);
        throw 1;
    }

    // Calculate the total number of grid cells.
    ntot = itot*jtot*ktot;

    // Calculate the grid dimensions per process.
    imax = itot / master->npx;
    jmax = jtot / master->npy;
    kmax = ktot;
    nmax = imax*jmax*kmax;

    // Calculate the block sizes for the transposes.
    iblock = itot / master->npy;
    jblock = jtot / master->npx;
    kblock = ktot / master->npx;

    // Calculate the grid dimensions including ghost cells.
    icells  = (imax+2*igc);
    jcells  = (jmax+2*jgc);
    ijcells = icells*jcells;
    kcells  = (kmax+2*kgc);
    ncells  = (imax+2*igc)*(jmax+2*jgc)*(kmax+2*kgc);

    // Calculate the starting and ending points for loops over the grid.
    istart = igc;
    jstart = jgc;
    kstart = kgc;

    iend   = imax + igc;
    jend   = jmax + jgc;
    kend   = kmax + kgc;

    check_ghost_cells();

    // allocate all arrays
    x     = new double[imax+2*igc];
    xh    = new double[imax+2*igc];
    y     = new double[jmax+2*jgc];
    yh    = new double[jmax+2*jgc];
    z     = new double[kmax+2*kgc];
    zh    = new double[kmax+2*kgc];
    dz    = new double[kmax+2*kgc];
    dzh   = new double[kmax+2*kgc];
    dzi   = new double[kmax+2*kgc];
    dzhi  = new double[kmax+2*kgc];
    dzi4  = new double[kmax+2*kgc];
    dzhi4 = new double[kmax+2*kgc];

    // allocate the data for the fourier transforms
    fftini  = fftw_alloc_real(itot*jmax);
    fftouti = fftw_alloc_real(itot*jmax);
    fftinj  = fftw_alloc_real(jtot*iblock);
    fftoutj = fftw_alloc_real(jtot*iblock);

    // initialize the communication functions
    init_mpi();
}

/**
 * This function initializes the fields containing the grid dimensions based
 * on the profiles in the input file.
 * @param inputin Pointer to the input class.
 */
void Grid::create(Input *inputin)
{
    // get the grid coordinates from the input
    if (inputin->get_prof(&z[kstart], "z", kmax))
        throw 1;

    if (z[kend-1] > zsize)
    {
        master->print_error("Highest grid point is above prescribed zsize\n");
        throw 1;
    }

    // calculate the grid
    calculate();
}

/**
 * This function calculates the scalars and arrays that contain the information
 * on the grid spacing.
 */
void Grid::calculate()
{
    int i,j,k;

    // calculate the grid spacing
    dx  = xsize / itot;
    dy  = ysize / jtot;
    dxi = 1./dx;
    dyi = 1./dy;

    // calculate the offset per process to get the true x- and y-coordinate
    double xoff = master->mpicoordx * xsize / master->npx;
    double yoff = master->mpicoordy * ysize / master->npy;

    // calculate the x and y coordinates
    for (i=0; i<icells; ++i)
    {
        x [i] = 0.5*dx + (i-igc)*dx + xoff;
        xh[i] = (i-igc)*dx + xoff;
    }

    for (j=0; j<jcells; ++j)
    {
        y [j] = 0.5*dy + (j-jgc)*dy + yoff;
        yh[j] = (j-jgc)*dy + yoff;
    }

    // the calculation of ghost cells and flux levels has to go according to numerical scheme
    if (swspatialorder == "2")
    {
        z[kstart-1] = -z[kstart];
        z[kend]     = 2.*zsize - z[kend-1];

        for (k=kstart+1; k<kend; ++k)
            zh[k] = 0.5*(z[k-1]+z[k]);
        zh[kstart] = 0.;
        zh[kend]   = zsize;

        // calculate the half levels according to the numerical scheme
        // compute the height of the grid cells
        for (k=1; k<kcells; ++k)
        {
            dzh [k] = z[k] - z[k-1];
            dzhi[k] = 1./dzh[k];
        }
        dzh [kstart-1] = dzh [kstart+1];
        dzhi[kstart-1] = dzhi[kstart+1];

        // compute the height of the grid cells
        for (k=1; k<kcells-1; ++k)
        {
            dz [k] = zh[k+1] - zh[k];
            dzi[k] = 1./dz[k];
        }
        dz [kstart-1] = dz [kstart];
        dzi[kstart-1] = dzi[kstart];
        dz [kend]     = dz [kend-1];
        dzi[kend]     = dzi[kend-1];

        // do not calculate 4th order gradients for 2nd order
    }

    if (swspatialorder == "4")
    {
        using namespace Finite_difference::O4;

        // calculate the height of the ghost cell
        z[kstart-1] = -2.*z[kstart] + (1./3.)*z[kstart+1];
        z[kstart-2] = -9.*z[kstart] +      2.*z[kstart+1];

        z[kend  ] = (8./3.)*zsize - 2.*z[kend-1] + (1./3.)*z[kend-2];
        z[kend+1] =      8.*zsize - 9.*z[kend-1] +      2.*z[kend-2];

        // Initialize the non-used values at a large value
        z[kstart-3] = Constants::dhuge;
        z[kend+2  ] = Constants::dhuge;

        zh[kstart] = 0.;
        for (k=kstart+1; k<kend; ++k)
            zh[k] = ci0*z[k-2] + ci1*z[k-1] + ci2*z[k] + ci3*z[k+1];
        zh[kend] = zsize;

        zh[kstart-1] = bi0*z[kstart-2] + bi1*z[kstart-1] + bi2*z[kstart] + bi3*z[kstart+1];
        zh[kend+1]   = ti0*z[kend-2  ] + ti1*z[kend-1  ] + ti2*z[kend  ] + ti3*z[kend+1  ];

        // calculate the half levels according to the numerical scheme
        // compute the height of the grid cells
        for (k=1; k<kcells; ++k)
        {
            dzh [k] = z[k] - z[k-1];
            dzhi[k] = 1./dzh[k];
        }
        dzh [kstart-3] = dzh [kstart+3];
        dzhi[kstart-3] = dzhi[kstart+3];

        // compute the height of the grid cells
        for (k=1; k<kcells-1; ++k)
        {
            dz [k] = zh[k+1] - zh[k];
            dzi[k] = 1./dz[k];
        }
        dz [kstart-3] = dz [kstart+2];
        dzi[kstart-3] = dzi[kstart+2];
        dz [kend+2] = dz [kend-3];
        dzi[kend+2] = dzi[kend-3];

        // calculate the fourth order gradients
        for (k=kstart; k<kend; ++k)
        {
            dzi4 [k] = 1./(cg0*zh[k-1] + cg1*zh[k  ] + cg2*zh[k+1] + cg3*zh[k+2]);
            dzhi4[k] = 1./(cg0*z [k-2] + cg1*z [k-1] + cg2*z [k  ] + cg3*z [k+1]);
        }
        dzhi4[kend  ] = 1./(cg0*z[kend-2] + cg1*z[kend-1] + cg2*z[kend] + cg3*z[kend+1]);

        // bc's
        dzi4 [kstart-1] = 1./(bg0*zh[kstart-1] + bg1*zh[kstart  ] + bg2*zh[kstart+1] + bg3*zh[kstart+2]);
        dzhi4[kstart-1] = 1./(bg0*z [kstart-2] + bg1*z [kstart-1] + bg2*z [kstart  ] + bg3*z [kstart+1]);

        dzi4 [kend  ] = 1./(tg0*zh[kend-2] + tg1*zh[kend-1] + tg2*zh[kend] + tg3*zh[kend+1]);
        dzhi4[kend+1] = 1./(tg0*z [kend-2] + tg1*z [kend-1] + tg2*z [kend] + tg3*z [kend+1]);

        // Define gradients at the boundary for the divgrad calculations.
        dzhi4bot = 1./(bg0*z[kstart-1] + bg1*z[kstart] + bg2*z[kstart+1] + bg3*z[kstart+2]);
        dzhi4top = 1./(tg0*z[kend-3  ] + tg1*z[kend-2] + tg2*z[kend-1  ] + tg3*z[kend    ]);

        // Initialize the unused values at a huge value to allow for easier error tracing.
        dzi4[kstart-2] = Constants::dhuge;
        dzi4[kstart-3] = Constants::dhuge;
        dzi4[kend+1  ] = Constants::dhuge;
        dzi4[kend+2  ] = Constants::dhuge;
    }
}

/**
 * This function checks whether the number of ghost cells does not exceed the slice thickness.
 */
void Grid::check_ghost_cells()
{
    // Check whether the size per patch is larger than number of ghost cells for 3D runs.
    if (imax < igc)
    {
	    master->print_error("Patch size in x-dir (%d) is smaller than the number of ghost cells (%d).\n",(iend-istart),igc);
	    master->print_error("Either increase itot or decrease npx.\n");
        throw 1;
    }

    // Check the jtot > 1 condition, to still allow for 2d runs.
    if (jtot > 1 && jmax < jgc)
    {
	    master->print_error("Patch size in y-dir (%d) is smaller than the number of ghost cells (%d).\n",(jend-jstart),jgc);
	    master->print_error("Either increase jtot or decrease npy.\n");
        throw 1;
    }
}

/**
 * This function increases the number of ghost cells in case necessary.
 * @param igc Ghost cells in the x-direction.
 * @param jgc Ghost cells in the y-direction.
 * @param kgc Ghost cells in the z-direction.
 */
void Grid::set_minimum_ghost_cells(const int igcin, const int jgcin, const int kgcin)
{
    igc = std::max(igc, igcin);
    jgc = std::max(jgc, jgcin);
    kgc = std::max(kgc, kgcin);

    // BvS: this doesn't work; imax is undefined if this routine is called from a class constructor
    // Removed it since this check is anyhow always performed from the init() of grid (after defining imax)
    //check_ghost_cells();
}

/**
 * This function does a second order horizontal interpolation in the x-direction
 * to the selected location on the grid.
 * @param out Pointer to the output field.
 * @param in Pointer to the input field.
 * @param locx Integer containing the location of the input field,
 * where a value of 1 refers to the flux level.
 */
void Grid::interpolate_2nd(double* const restrict out, const double* const restrict in, const int locin[3], const int locout[3])
{
    const int ii = 1;
    const int jj = icells;
    const int kk = ijcells;

    const int iih = (locin[0]-locout[0])*ii;
    const int jjh = (locin[1]-locout[1])*jj;

    // interpolate the field
    // \TODO add the vertical component
    for (int k=0; k<kcells; ++k)
        for (int j=jstart; j<jend; ++j)
#pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                out[ijk] = 0.5*(0.5*in[ijk    ] + 0.5*in[ijk+iih    ])
                         + 0.5*(0.5*in[ijk+jjh] + 0.5*in[ijk+iih+jjh]);
            }
}

/**
 * This function does a fourth order horizontal interpolation in the x-direction
 * to the selected location on the grid.
 * @param out Pointer to the output field.
 * @param in Pointer to the input field.
 * @param locx Integer containing the location of the input field,
 * where a value of 1 refers to the flux level.
 */
void Grid::interpolate_4th(double* restrict out, double* restrict in, const int locin[3], const int locout[3])
{
    using namespace Finite_difference::O4;

    // interpolation function, locx = 1 indicates that the reference is at the half level
    const int ii = 1;
    const int jj = icells;
    const int kk = ijcells;

    // a shift to the left gives minus 1 a shift to the right +1
    const int iih1 = 1*(locin[0]-locout[0])*ii;
    const int iih2 = 2*(locin[0]-locout[0])*ii;
    const int jjh1 = 1*(locin[1]-locout[1])*jj;
    const int jjh2 = 2*(locin[1]-locout[1])*jj;

    // \TODO add the vertical component
    for (int k=0; k<kcells; ++k)
        for (int j=jstart; j<jend; ++j)
#pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                out[ijk] = ci0*(ci0*in[ijk-iih1-jjh1] + ci1*in[ijk-jjh1] + ci2*in[ijk+iih1-jjh1] + ci3*in[ijk+iih2-jjh1])
                         + ci1*(ci0*in[ijk-iih1     ] + ci1*in[ijk     ] + ci2*in[ijk+iih1     ] + ci3*in[ijk+iih2     ])
                         + ci2*(ci0*in[ijk-iih1+jjh1] + ci1*in[ijk+jjh1] + ci2*in[ijk+iih1+jjh1] + ci3*in[ijk+iih2+jjh1])
                         + ci3*(ci0*in[ijk-iih1+jjh2] + ci1*in[ijk+jjh2] + ci2*in[ijk+iih1+jjh2] + ci3*in[ijk+iih2+jjh2]);
            }
}

void Grid::calc_mean(double* restrict prof, const double* restrict data, const int krange)
{
    const int jj = icells;
    const int kk = ijcells;

    for (int k=0; k<krange; ++k)
    {
        prof[k] = 0.;
        for (int j=jstart; j<jend; ++j)
#pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk  = i + j*jj + k*kk;
                prof[k] += data[ijk];
            }
    }

    master->sum(prof, krange);

    const double n = itot*jtot;

    for (int k=0; k<krange; ++k)
        prof[k] /= n;
}
