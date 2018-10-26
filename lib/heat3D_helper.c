/*
 * heat3D_helper.c
 *
 *  Created on: Sep 26, 2018
 *      Author: Derek W. Harrison
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../inc/memory_functions.h"
#include "../inc/mappings.h"
#include "../inc/user_types.h"


/*-----------------------------------------------------------------------------------------------*/
void vector_addition(double *v1,
                     double f_v1,
                     double *v2,
                     double f_v2,
                     int size_vec,
                     double *vr)
{
    /*
     * Addition of vectors:
     *
     * vr = f_v1*v1 + f_v2*v2
     *
     * input    v1
     * input    factor_v1
     * input    v2
     * input    factor_v2
     * input    size_vector
     * output   output_vector
     */

    int j;

    for (j = 1; j <= size_vec; j++)
        vr[j] = f_v1*v1[j] + f_v2*v2[j];

}


/*-----------------------------------------------------------------------------------------------*/
void dot_product(double *v1,
                 double *v2,
                 int size_vec,
                 double *vr)
{
    /*
     * Calculating the dot product:
     *
     * v_r = v_1'*v_2
     *
     * input    v_1
     * input    v_2
     * input    size_vec
     * output   v_r
     */

    int j;

    double dummy = 0.0;
    for (j = 1; j <= size_vec; j++)
        dummy = dummy + v1[j]*v2[j];

    *vr = dummy;

}


/*-----------------------------------------------------------------------------------------------*/
void mat_vec_mult(grid_size_t grid_size,
                  double** A,
                  double* p,
                  double* Ap)
{
    /*
     * Calculating the matrix product A*p where A is the vectorized coefficient matrix:
     *
     * Ap = A*p
     *
     * input    grid_size
     * input    A
     * input    p
     * output   Ap
     */

    int j;
    int nt;

    nt = grid_size.nx*grid_size.ny*grid_size.nz;

    Ap[1] = A[1][4]*p[1] + A[2][3]*p[2] + A[1+grid_size.nx][2]*p[1+grid_size.nx] +
            A[1+grid_size.nx*grid_size.ny][1]*p[1+grid_size.nx*grid_size.ny];

    for (j = 2; j <= grid_size.nx; j++)
        Ap[j] = A[j][3]*p[j-1] + A[j][4]*p[j] + A[j+1][3]*p[j+1] +
        A[j+grid_size.nx][2]*p[j+grid_size.nx] +
        A[j+grid_size.nx*grid_size.ny][1]*p[j+grid_size.nx*grid_size.ny];

    for (j= grid_size.nx + 1; j <= grid_size.nx*grid_size.ny; j++)
        Ap[j] = A[j][2]*p[j-grid_size.nx] + A[j][3]*p[j-1] + A[j][4]*p[j] + A[j+1][3]*p[j+1] +
        A[j+grid_size.nx][2]*p[j+grid_size.nx] +
        A[j+grid_size.nx*grid_size.ny][1]*p[j+grid_size.nx*grid_size.ny];

    for (j = grid_size.nx*grid_size.ny + 1; j <= nt - grid_size.nx*grid_size.ny; j++)
        Ap[j] = A[j][1]*p[j-grid_size.nx*grid_size.ny] + A[j][2]*p[j-grid_size.nx] +
        A[j][3]*p[j-1] + A[j][4]*p[j] + A[j+1][3]*p[j+1] + A[j+grid_size.nx][2]*p[j+grid_size.nx] +
        A[j+grid_size.nx*grid_size.ny][1]*p[j+grid_size.nx*grid_size.ny];

    for (j=nt - grid_size.nx*grid_size.ny + 1; j <= nt - grid_size.nx; j++)
        Ap[j] = A[j][1]*p[j-grid_size.nx*grid_size.ny] + A[j][2]*p[j-grid_size.nx] +
        A[j][3]*p[j-1] + A[j][4]*p[j] + A[j+1][3]*p[j+1] + A[j+grid_size.nx][2]*p[j+grid_size.nx];

    for (j=nt - grid_size.nx + 1; j <= nt - 1; j++)
        Ap[j] = A[j][1]*p[j-grid_size.nx*grid_size.ny] + A[j][2]*p[j-grid_size.nx] +
        A[j][3]*p[j-1] + A[j][4]*p[j] + A[j+1][3]*p[j+1];

    Ap[nt] = A[nt][1]*p[nt-grid_size.nx*grid_size.ny] + A[nt][2]*p[nt-grid_size.nx] + A
            [nt][3]*p[nt-1] + A[nt][4]*p[nt];

}


/*-----------------------------------------------------------------------------------------------*/
void generate_coefficient_matrix(domain_size_t domain_size,
                                 grid_size_t grid_size,
                                 boundary_conditions_t boundary_conditions,
                                 time_data_t time_data,
                                 physical_paramaters_t physical_parameters,
                                 double (*source)(double x,double y,double z,double t),
                                 grid_coordinates_t* grid_coordinates,
                                 kershaw_algorithm_data_t* kershaw_data)
{
    /*
     * Generated the vectorized coefficient matrix for the poisson solver
     *
     * input        domain_size
     * input        grid_size
     * input        boundary_conditions
     * input        time_data
     * input        physical_parameters
     * input        source
     * input        grid_coordinates
     * input/output kershaw_data
     */

    double deltax, deltay, deltaz;
    double b1, b2, b3;
    double t, K;
    double ***X, ***Y, ***Z;
    double *xo, *r, **A;
    int nn;
    int i, j, k;
    double *wb_switch, *eb_switch, *sb_switch;
    double *nb_switch, *bb_switch, *tb_switch;
    boundary_type_t wb_type, eb_type, sb_type;
    boundary_type_t nb_type, bb_type, tb_type;
    boundary_funcs_t boundary_funcs;


    /* Setting boundries and boundary types */
    boundary_funcs = boundary_conditions.boundary_funcs;

    wb_type = boundary_conditions.boundary_type_faces.west_boundary;
    eb_type = boundary_conditions.boundary_type_faces.east_boundary;
    sb_type = boundary_conditions.boundary_type_faces.south_boundary;
    nb_type = boundary_conditions.boundary_type_faces.north_boundary;
    bb_type = boundary_conditions.boundary_type_faces.bottom_boundary;
    tb_type = boundary_conditions.boundary_type_faces.top_boundary;


    /* Initializing parameters */
    t  = time_data.t;
    K  = physical_parameters.rho*physical_parameters.Cp/time_data.dt;

    X = grid_coordinates->X;
    Y = grid_coordinates->Y;
    Z = grid_coordinates->Z;

    deltax = domain_size.Lx/grid_size.nx;
    deltay = domain_size.Ly/grid_size.ny;
    deltaz = domain_size.Lz/grid_size.nz;

    b1 = -physical_parameters.conductivity.gammax/(deltax*deltax);
    b2 = -physical_parameters.conductivity.gammay/(deltay*deltay);
    b3 = -physical_parameters.conductivity.gammaz/(deltaz*deltaz);

    xo = kershaw_data->x;
    r  = kershaw_data->r;
    A  = kershaw_data->A;


    /* Setting switches */
    wb_switch = boundary_switch_mapper(deltax * 2 * b1);
    eb_switch = boundary_switch_mapper(deltax * 2 * b1);
    sb_switch = boundary_switch_mapper(deltay * 2 * b2);
    nb_switch = boundary_switch_mapper(deltay * 2 * b2);
    bb_switch = boundary_switch_mapper(deltaz * 2 * b3);
    tb_switch = boundary_switch_mapper(deltaz * 2 * b3);


    /* Generating vectorized coefficient matrix */
    //Generating central coefficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = -2*b1 - 2*b2 - 2*b3 + K;
        r[nn] = source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + b2*xo[nn-grid_size.nx] + b1*xo[nn-1]
                + A[nn][4]*xo[nn] + b1*xo[nn+1] + b2*xo[nn+grid_size.nx]
                + b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating corner 1 coefficients
    for (i = 1; i <= 1; i++)
    for (j = 1; j <= 1; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = 0;
        A[nn][3] = 0;
        A[nn][4] = -(b1 + 2*b1*wb_type) - (b2 + 2*b2*sb_type) - (b3 + 2*b3*bb_type) + K;
        r[nn] = - 2*b1*wb_switch[wb_type]*boundary_funcs.boundary_west(Y[i][j][k], Z[i][j][k], t)
                - 2*b2*sb_switch[sb_type]*boundary_funcs.boundary_south(X[i][j][k], Z[i][j][k], t)
                - 2*b3*bb_switch[bb_type]*boundary_funcs.boundary_bottom(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (A[nn][4]*xo[nn] + b1*xo[nn+1] + b2*xo[nn+grid_size.nx] +
                   b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating side a cofficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = 1; j <= 1; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = 0;
        A[nn][3] = b1;
        A[nn][4] = -2*b1 - (b2 + 2*b2*sb_type) - (b3 + 2*b3*bb_type) + K;
        r[nn] = - 2*b2*sb_switch[sb_type]*boundary_funcs.boundary_south(X[i][j][k], Z[i][j][k], t)
                - 2*b3*bb_switch[bb_type]*boundary_funcs.boundary_bottom(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b1*xo[nn-1] + A[nn][4]*xo[nn] + b1*xo[nn+1] + b2*xo[nn+grid_size.nx] +
                   b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating corner 2 coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = 1; j <= 1; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = 0;
        A[nn][3] = b1;
        A[nn][4] = -(b1 + 2*b1*eb_type) - (b2 + 2*b2*sb_type) - (b3 + 2*b3*bb_type) + K;
        r[nn] = - 2*b1*eb_switch[eb_type]*boundary_funcs.boundary_east(Y[i][j][k], Z[i][j][k], t)
                - 2*b2*sb_switch[sb_type]*boundary_funcs.boundary_south(X[i][j][k], Z[i][j][k], t)
                - 2*b3*bb_switch[bb_type]*boundary_funcs.boundary_bottom(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b1*xo[nn-1] + A[nn][4]*xo[nn] + b2*xo[nn+grid_size.nx] +
                   b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating side d coefficients and source terms
    for (i = 1; i <= 1; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = b2;
        A[nn][3] = 0;
        A[nn][4] = -(b1 + 2*b1*wb_type) - 2*b2 - (b3 + 2*b3*bb_type) + K;
        r[nn] = - 2*b1*wb_switch[wb_type]*boundary_funcs.boundary_west(Y[i][j][k], Z[i][j][k], t)
                - 2*b3*bb_switch[bb_type]*boundary_funcs.boundary_bottom(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b2*xo[nn-grid_size.nx] + A[nn][4]*xo[nn] + b1*xo[nn+1] + b2*xo[nn+grid_size.nx] +
                   b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating face E coefficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = -2*b1 - 2*b2 - (b3 + 2*b3*bb_type) + K;
        r[nn] = - 2*b3*bb_switch[bb_type]*boundary_funcs.boundary_bottom(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b2*xo[nn-grid_size.nx] + b1*xo[nn-1] + A[nn][4]*xo[nn] + b1*xo[nn+1] +
                   b2*xo[nn+grid_size.nx] + b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating side b coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = -(b1 + 2*b1*eb_type) - 2*b2 - (b3 + 2*b3*bb_type) + K;
        r[nn] = - 2*b1*eb_switch[eb_type]*boundary_funcs.boundary_east(Y[i][j][k], Z[i][j][k], t)
                - 2*b3*bb_switch[bb_type]*boundary_funcs.boundary_bottom(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b2*xo[nn-grid_size.nx] + b1*xo[nn-1] + A[nn][4]*xo[nn] + b2*xo[nn+grid_size.nx] +
                   b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating corner 4 coefficients and source terms
    for (i = 1; i <= 1; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = b2;
        A[nn][3] = 0;
        A[nn][4] = -(b1 + 2*b1*wb_type) - (b2 + 2*b2*nb_type) - (b3 + 2*b3*bb_type) + K;
        r[nn] = - 2*b1*wb_switch[wb_type]*boundary_funcs.boundary_west(Y[i][j][k], Z[i][j][k], t)
                - 2*b2*nb_switch[nb_type]*boundary_funcs.boundary_north(X[i][j][k], Z[i][j][k], t)
                - 2*b3*bb_switch[bb_type]*boundary_funcs.boundary_bottom(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b2*xo[nn-grid_size.nx] + A[nn][4]*xo[nn] + b1*xo[nn+1] +
                   b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating side c coefficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = -2*b1 - (b2 + 2*b2*nb_type) - (b3 + 2*b3*bb_type) + K;
        r[nn] = - 2*b2*nb_switch[nb_type]*boundary_funcs.boundary_north(X[i][j][k], Z[i][j][k], t)
                - 2*b3*bb_switch[bb_type]*boundary_funcs.boundary_bottom(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b2*xo[nn-grid_size.nx] + b1*xo[nn-1] + A[nn][4]*xo[nn] + b1*xo[nn+1] +
                   b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating corner 3 coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = 1; k <= 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = 0;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = -(b1 + 2*b1*eb_type) - (b2 + 2*b2*nb_type) - (b3 + 2*b3*bb_type) + K;
        r[nn] = - 2*b1*eb_switch[eb_type]*boundary_funcs.boundary_east(Y[i][j][k], Z[i][j][k], t)
                - 2*b2*nb_switch[nb_type]*boundary_funcs.boundary_north(X[i][j][k], Z[i][j][k], t)
                - 2*b3*bb_switch[bb_type]*boundary_funcs.boundary_bottom(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b2*xo[nn-grid_size.nx] + b1*xo[nn-1] + A[nn][4]*xo[nn] +
                   b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating side e coefficients and source terms
    for (i = 1; i <= 1; i++)
    for (j = 1; j <= 1; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = 0;
        A[nn][3] = 0;
        A[nn][4] = -(b1 + 2*b1*wb_type) - (b2 + 2*b2*sb_type) - 2*b3 + K;
        r[nn] = - 2*b1*wb_switch[wb_type]*boundary_funcs.boundary_west(Y[i][j][k], Z[i][j][k], t)
                - 2*b2*sb_switch[sb_type]*boundary_funcs.boundary_south(X[i][j][k], Z[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + A[nn][4]*xo[nn] + b1*xo[nn+1] +
                   b2*xo[nn+grid_size.nx] + b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating face A coefficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = 1; j <= 1; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = 0;
        A[nn][3] = b1;
        A[nn][4] = -2*b1 - (b2 + 2*b2*sb_type) - 2*b3 + K;
        r[nn] = - 2*b2*sb_switch[sb_type]*boundary_funcs.boundary_south(X[i][j][k], Z[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + b1*xo[nn-1] + A[nn][4]*xo[nn] + b1*xo[nn+1] +
                   b2*xo[nn+grid_size.nx] + b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating side f coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = 1; j <= 1; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = 0;
        A[nn][3] = b1;
        A[nn][4] = -(b1 + 2*b1*eb_type) - (b2 + 2*b2*sb_type) - 2*b3 + K;
        r[nn] = - 2*b1*eb_switch[eb_type]*boundary_funcs.boundary_east(Y[i][j][k], Z[i][j][k], t)
                - 2*b2*sb_switch[sb_type]*boundary_funcs.boundary_south(X[i][j][k], Z[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + b1*xo[nn-1] + A[nn][4]*xo[nn] +
                   b2*xo[nn+grid_size.nx] + b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating face D coefficients and source terms
    for (i = 1; i <= 1; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = 0;
        A[nn][4] = -(b1 + 2*b1*wb_type) - 2*b2 - 2*b3 + K;
        r[nn] = - 2*b1*wb_switch[wb_type]*boundary_funcs.boundary_west(Y[i][j][k], Z[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + b2*xo[nn-grid_size.nx] + A[nn][4]*xo[nn] +
                   b1*xo[nn+1] + b2*xo[nn+grid_size.nx] + b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating face B coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = -(b1 + 2*b1*eb_type) - 2*b2 - 2*b3 + K;
        r[nn] = - 2*b1*eb_switch[eb_type]*boundary_funcs.boundary_east(Y[i][j][k], Z[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + b2*xo[nn-grid_size.nx] + b1*xo[nn-1] +
                   A[nn][4]*xo[nn] + b2*xo[nn+grid_size.nx] + b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating side h coefficients and source terms
    for (i = 1; i <= 1; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = 0;
        A[nn][4] = -(b1 + 2*b1*wb_type) - (b2 + 2*b2*nb_type) - 2*b3 + K;
        r[nn] = - 2*b1*wb_switch[wb_type]*boundary_funcs.boundary_west(Y[i][j][k], Z[i][j][k], t)
                - 2*b2*nb_switch[nb_type]*boundary_funcs.boundary_north(X[i][j][k], Z[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + b2*xo[nn-grid_size.nx] + A[nn][4]*xo[nn] +
                   b1*xo[nn+1] + b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating face C coefficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = -2*b1 - (b2 + 2*b2*nb_type) - 2*b3 + K;
        r[nn] = - 2*b2*nb_switch[nb_type]*boundary_funcs.boundary_north(X[i][j][k], Z[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + b2*xo[nn-grid_size.nx] + b1*xo[nn-1] +
                   A[nn][4]*xo[nn] + b1*xo[nn+1] + b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating side g coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = 2; k <= grid_size.nz - 1; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = -(b1 + 2*b1*eb_type) - (b2 + 2*b2*nb_type) - 2*b3 + K;
        r[nn] = - 2*b1*eb_switch[eb_type]*boundary_funcs.boundary_east(Y[i][j][k], Z[i][j][k], t)
                - 2*b2*nb_switch[nb_type]*boundary_funcs.boundary_north(X[i][j][k], Z[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + b2*xo[nn-grid_size.nx] + b1*xo[nn-1] +
                   A[nn][4]*xo[nn] + b3*xo[nn+grid_size.nx*grid_size.ny]);
    }

    //Generating corner 5 coefficients and source terms
    for (i = 1; i <= 1; i++)
    for (j = 1; j <= 1; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = 0;
        A[nn][3] = 0;
        A[nn][4] = -(b1 + 2*b1*wb_type) - (b2 + 2*b2*sb_type) - (b3 + 2*b3*tb_type) + K;
        r[nn] = - 2*b1*wb_switch[wb_type]*boundary_funcs.boundary_west(Y[i][j][k], Z[i][j][k], t)
                - 2*b2*sb_switch[sb_type]*boundary_funcs.boundary_south(X[i][j][k], Z[i][j][k], t)
                - 2*b3*tb_switch[tb_type]*boundary_funcs.boundary_top(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + A[nn][4]*xo[nn] + b1*xo[nn+1] +
                   b2*xo[nn+grid_size.nx]);
    }

    //Generating side i coefficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = 1; j <= 1; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = 0;
        A[nn][3] = b1;
        A[nn][4] = -2*b1 - (b2 + 2*b2*sb_type) - (b3 + 2*b3*tb_type) + K;
        r[nn] = - 2*b2*sb_switch[sb_type]*boundary_funcs.boundary_south(X[i][j][k], Z[i][j][k], t)
                - 2*b3*tb_switch[tb_type]*boundary_funcs.boundary_top(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + b1*xo[nn-1] + A[nn][4]*xo[nn] +
                   b1*xo[nn+1] + b2*xo[nn+grid_size.nx]);
    }

    //Generating corner 6 coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = 1; j <= 1; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = 0;
        A[nn][3] = b1;
        A[nn][4] = -(b1 + 2*b1*eb_type) - (b2 + 2*b2*sb_type) - (b3 + 2*b3*tb_type) + K;
        r[nn] = - 2*b1*eb_switch[eb_type]*boundary_funcs.boundary_east(Y[i][j][k], Z[i][j][k], t)
                - 2*b2*sb_switch[sb_type]*boundary_funcs.boundary_south(X[i][j][k], Z[i][j][k], t)
                - 2*b3*tb_switch[tb_type]*boundary_funcs.boundary_top(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + b1*xo[nn-1] + A[nn][4]*xo[nn] +
                   b2*xo[nn+grid_size.nx]);
    }

    //Generating side l coefficients and source terms
    for (i = 1; i <= 1; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = 0;
        A[nn][4] = -(b1 + 2*b1*wb_type) - 2*b2 - (b3 + 2*b3*tb_type) + K;
        r[nn] = - 2*b1*wb_switch[wb_type]*boundary_funcs.boundary_west(Y[i][j][k], Z[i][j][k], t)
                - 2*b3*tb_switch[tb_type]*boundary_funcs.boundary_top(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + b2*xo[nn-grid_size.nx] + A[nn][4]*xo[nn] +
                   b1*xo[nn+1] + b2*xo[nn+grid_size.nx]);
    }

    //Generating face F coefficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = -2*b1 - 2*b2 - (b3 + 2*b3*tb_type) + K;
        r[nn] = - 2*b3*tb_switch[tb_type]*boundary_funcs.boundary_top(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + b2*xo[nn-grid_size.nx] + b1*xo[nn-1] +
                   A[nn][4]*xo[nn] + b1*xo[nn+1] + b2*xo[nn+grid_size.nx]);
    }

    //Generating side j coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = 2; j <= grid_size.ny - 1; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = -(b1 + 2*b1*eb_type) - 2*b2 - (b3 + 2*b3*tb_type) + K;
        r[nn] = - 2*b1*eb_switch[eb_type]*boundary_funcs.boundary_east(Y[i][j][k], Z[i][j][k], t)
                - 2*b3*tb_switch[tb_type]*boundary_funcs.boundary_top(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + b2*xo[nn-grid_size.nx] + b1*xo[nn-1] +
                   A[nn][4]*xo[nn] + b2*xo[nn+grid_size.nx]);
    }

    //Generating corner 8 coefficients and source terms
    for (i = 1; i <= 1; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = 0;
        A[nn][4] = -(b1 + 2*b1*wb_type) - (b2 + 2*b2*nb_type) - (b3 + 2*b3*tb_type) + K;
        r[nn] = - 2*b1*wb_switch[wb_type]*boundary_funcs.boundary_west(Y[i][j][k], Z[i][j][k], t)
                - 2*b2*nb_switch[nb_type]*boundary_funcs.boundary_north(X[i][j][k], Z[i][j][k], t)
                - 2*b3*tb_switch[tb_type]*boundary_funcs.boundary_top(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + b2*xo[nn-grid_size.nx] + A[nn][4]*xo[nn] +
                   b1*xo[nn+1]);
    }

    //Generating side k coefficients and source terms
    for (i = 2; i <= grid_size.nx - 1; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = -2*b1 - (b2 + 2*b2*nb_type) - (b3 + 2*b3*tb_type) + K;
        r[nn] = - 2*b2*nb_switch[nb_type]*boundary_funcs.boundary_north(X[i][j][k], Z[i][j][k], t)
                - 2*b3*tb_switch[tb_type]*boundary_funcs.boundary_top(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + b2*xo[nn-grid_size.nx] + b1*xo[nn-1] +
                   A[nn][4]*xo[nn] + b1*xo[nn+1]);
    }

    //Generating corner 7 coefficients and source terms
    for (i = grid_size.nx; i <= grid_size.nx; i++)
    for (j = grid_size.ny; j <= grid_size.ny; j++)
    for (k = grid_size.nz; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        A[nn][1] = b3;
        A[nn][2] = b2;
        A[nn][3] = b1;
        A[nn][4] = -(b1 + 2*b1*eb_type) - (b2 + 2*b2*nb_type) - (b3 + 2*b3*tb_type) + K;
        r[nn] = - 2*b1*eb_switch[eb_type]*boundary_funcs.boundary_east(Y[i][j][k], Z[i][j][k], t)
                - 2*b2*nb_switch[nb_type]*boundary_funcs.boundary_north(X[i][j][k], Z[i][j][k], t)
                - 2*b3*tb_switch[tb_type]*boundary_funcs.boundary_top(X[i][j][k], Y[i][j][k], t)
                + source(X[i][j][k], Y[i][j][k], Z[i][j][k], t) + K*xo[nn]
                - (b3*xo[nn-grid_size.nx*grid_size.ny] + b2*xo[nn-grid_size.nx] + b1*xo[nn-1] +
                   A[nn][4]*xo[nn]);
    }


    /* Deallocate memory */
    free(wb_switch);
    free(eb_switch);
    free(sb_switch);
    free(nb_switch);
    free(bb_switch);
    free(tb_switch);

}


/*-----------------------------------------------------------------------------------------------*/
void Ly_solver(grid_size_t grid_size,
               double** L,
               double* r,
               double* y)
{
    /*
     * Solves the linear system Ly = r where L is lower triangular
     * and L trivially vectorized
     *
     * input    grid_size
     * input    L
     * input    r
     * output   y
     */

    int j;
    int nt;

    nt = grid_size.nx*grid_size.ny*grid_size.nz;

    y[1] = r[1]/L[1][4];
    for (j = 2; j <= grid_size.nx; j++)
        y[j] = (r[j]-L[j][3]*y[j-1])/L[j][4];

    for (j = grid_size.nx + 1; j <= grid_size.nx*grid_size.ny; j++)
        y[j] = (r[j] - L[j][2]*y[j-grid_size.nx] - L[j][3]*y[j-1])/L[j][4];

    for (j = grid_size.nx*grid_size.ny + 1; j <= nt; j++)
        y[j] = (r[j] - L[j][1]*y[j-grid_size.nx*grid_size.ny] - L[j][2]*y[j-grid_size.nx] -
                L[j][3]*y[j-1])/L[j][4];

}


/*-----------------------------------------------------------------------------------------------*/
void LTz_solver(grid_size_t grid_size,
                double** L,
                double* y,
                double* z)
{
    /*
     * Solves the linear system Ly = z where L is upper triangular
     * and L trivially vectorized
     *
     * input    grid_size
     * input    L
     * input    y
     * output   z
     */

    int j;
    int nt;

    nt = grid_size.nx*grid_size.ny*grid_size.nz;

    z[nt] = y[nt]/L[nt][4];
    for(j = nt - 1; j >= nt - grid_size.nx + 1; j--)
        z[j] = (y[j] - L[j+1][3]*z[j+1])/L[j][4];

    for(j = nt - grid_size.nx; j >= nt - grid_size.nx*grid_size.ny + 1; j--)
        z[j] = (y[j] - L[j+1][3]*z[j+1] - L[j+grid_size.nx][2]*z[j+grid_size.nx])/L[j][4];

    for(j = nt - grid_size.nx*grid_size.ny; j >= 1; j--)
        z[j] = (y[j] - L[j+1][3]*z[j+1] - L[j+grid_size.nx][2]*z[j+grid_size.nx] -
                L[j+grid_size.nx*grid_size.ny][1]*z[j+grid_size.nx*grid_size.ny])/L[j][4];

}


/*-----------------------------------------------------------------------------------------------*/
void incomplete_cholesky_factorization(grid_size_t grid_size,
                                       double** A,
                                       double** L)
{
    /*
     * Performs the incomplete cholesky factorization on vectorized matrix A
     *
     * input    grid_size.
     * input     A
     * output    L
     */

    int i,j;
    int nt;


    nt = grid_size.nx*grid_size.ny*grid_size.nz;

    L[1][1] = 0;
    L[1][2] = 0;
    L[1][3] = 0;
    L[1][4] = sqrt(A[1][4]);

    for (j = 2; j <= nt; j++)
    for (i = 1; i <= 4; i++)
        if (A[j][i] != 0)
        {
            if (j>=1+grid_size.nx*grid_size.ny)
                L[j][1]=A[j][1]/L[j-grid_size.nx*grid_size.ny][4];
            if (j>=1+grid_size.nx)
                L[j][2]=A[j][2]/L[j-grid_size.nx][4];
            L[j][3]=A[j][3]/L[j-1][4];
            L[j][4]=sqrt(A[j][4] - L[j][1]*L[j][1] - L[j][2]*L[j][2] - L[j][3]*L[j][3]);
        }

}


void initialize_temperature_field(grid_size_t grid_size,
                                  time_data_t time_data,
                                  double *temp_field)
{
    /*
     * Initialize temperature field
     *
     * input    grid_size.
     * input    time_data
     * output   x
     */

    int j, nt;

    nt = grid_size.nx*grid_size.ny*grid_size.nz;

    for (j = 1; j <= nt; j++)
        temp_field[j]  = time_data.Tinitial;

}


/*-----------------------------------------------------------------------------------------------*/
void initialize_time_data(time_data_t* time_data)
{
    /*
     * Initialize temperature field
     *
     * input/output   time_data
     */

    time_data->t = time_data->ti;
    time_data->current_timestep = 0;
    time_data->dt = (time_data->tf - time_data->ti)/time_data->timesteps;

}


/*-----------------------------------------------------------------------------------------------*/
void generate_grid_coordinates(domain_size_t domain_size,
                               grid_size_t grid_size,
                               grid_coordinates_t* grid_coordinates)
{
    /*
     * Generate grid coordinates
     *
     * input    domain_size.
     * input    grid_size
     * output   grid_coordinates
     */

    int i, j, k;
    double deltax, deltay, deltaz;

    deltax = domain_size.Lx/grid_size.nx;
    deltay = domain_size.Ly/grid_size.ny;
    deltaz = domain_size.Lz/grid_size.nz;

    /*Generating node coordinates*/
    for (i = 1; i <= grid_size.nx; i++)
    for (j = 1; j <= grid_size.ny; j++)
    for (k = 1; k <= grid_size.nz; k++)
    {
        grid_coordinates->X[i][j][k] = i*deltax-deltax/2;
        grid_coordinates->Y[i][j][k] = j*deltay-deltay/2;
        grid_coordinates->Z[i][j][k] = k*deltaz-deltaz/2;
    }

}


/*-----------------------------------------------------------------------------------------------*/
void preconditioning(grid_size_t grid_size,
                     kershaw_algorithm_data_t* kershaw_data)
{
    /*
     * Precondition coefficient matrix A via incomplete cholesky factorization
     *
     * input    grid_size.
     * input    A
     * input    r
     * output   y
     * output   L
     */

    incomplete_cholesky_factorization(grid_size,
                                      kershaw_data->A,
                                      kershaw_data->L);

}


/*-----------------------------------------------------------------------------------------------*/
void execute_kershaw_algorithm(grid_size_t grid_size,
                               kershaw_algorithm_data_t* kershaw_data)
{
    /*
     * Execute kershaw algorithm
     *
     * input        grid_size
     * input/output kershaw_data
     */

    double delold, delnew, pAp, error;
    double alpha, B;
    int nt, imax;

    /*Initializing parameters*/
    imax  = 5000;
    error = 1e-30;
    nt    = grid_size.nx*grid_size.ny*grid_size.nz;

    dot_product(kershaw_data->r,
                kershaw_data->r,
                nt,
                &(kershaw_data->epsilon));

    Ly_solver(grid_size,
              kershaw_data->L,
              kershaw_data->r,
              kershaw_data->y);

    LTz_solver(grid_size,
               kershaw_data->L,
               kershaw_data->y,
               kershaw_data->p);

    kershaw_data->iterations = 0;
    do
    {
        dot_product(kershaw_data->r,
                    kershaw_data->p,
                    nt,
                    &delold);

        mat_vec_mult(grid_size,
                     kershaw_data->A,
                     kershaw_data->p,
                     kershaw_data->Ap);

        dot_product(kershaw_data->p,
                    kershaw_data->Ap,
                    nt,
                    &pAp);

        alpha = delold/pAp;

        vector_addition(kershaw_data->x,
                        1.0,
                        kershaw_data->p,
                        alpha,
                        nt,
                        kershaw_data->x);

        vector_addition(kershaw_data->r,
                        1.0,
                        kershaw_data->Ap,
                        -alpha,
                        nt,
                        kershaw_data->r);

        Ly_solver(grid_size,
                  kershaw_data->L,
                  kershaw_data->r,
                  kershaw_data->y);

        LTz_solver(grid_size,
                   kershaw_data->L,
                   kershaw_data->y,
                   kershaw_data->z);

        dot_product(kershaw_data->r,
                    kershaw_data->z,
                    nt,
                    &delnew);

        B = delnew/delold;

        vector_addition(kershaw_data->z,
                        1.0,
                        kershaw_data->p,
                        B,
                        nt,
                        kershaw_data->p);

        dot_product(kershaw_data->r,
                    kershaw_data->r,
                    nt,
                    &(kershaw_data->epsilon));

        kershaw_data->epsilon = sqrt(kershaw_data->epsilon/nt);
        kershaw_data->iterations = kershaw_data->iterations + 1;

    }while (kershaw_data->iterations < imax && kershaw_data->epsilon > error);

}


/*-----------------------------------------------------------------------------------------------*/
void processing_results(grid_size_t grid_size,
                        kershaw_algorithm_data_t* kershaw_data,
                        double*** T)
{
    /*
     * Process results temperature field
     *
     * input    grid_size
     * input    kershaw_data
     * output   T
     */

    int i, j, k, nn;

    for (i = 1; i <= grid_size.nx ; i++)
    for (j = 1; j <= grid_size.ny ; j++)
    for (k = 1; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        T[i][j][k] = kershaw_data->x[nn];
    }

}
