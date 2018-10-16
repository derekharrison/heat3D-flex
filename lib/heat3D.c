/*
 * heat3D.c
 *
 *  Created on: Sep 25, 2018
 *      Author: Derek W. Harrison
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../inc/heat3D_helper.h"
#include "../inc/memory_functions.h"
#include "../inc/user_types.h"


/*-----------------------------------------------------------------------------------------------*/
void heat3D(domain_size_t domain_size,
            grid_size_t grid_size,
            boundary_conditions_t boundary_conditions,
            time_dep_input_t time_dep_input,
            gammas_t gammas,
            double (*source)(double x,double y,double z,double t),
            grid_coordinates_t* grid_coordinates,
            double ***T)
{
    /*
     * This function solves the 3D heat equation:
     *
     * gammax*d2T/dx2 + gammay*d2T/dy2 + gammaz*d2T/dz2 + q(x,y,z,t) = rho*Cp*dT/dt
     *
     * on a rectangluar grid. The code handles mixed boundary conditions, i.e.
     * both Neumann and Dirichlet boundary conditions.
     *
     * The equation is solved using a preconditioned conjugate gradient method with
     * incomplete cholesky factorization as preconditioner.
     *
     * Note: a node numbering scheme was used such that the first element of all
     * arrays start at 1 as opposed to zero. The numbering scheme used in this work is
     *
     * nn = i + (j-1)*nx + (k-1)*nx*ny
     *
     * where nn is the node number, nx the number of nodes used in the x direction,
     * ny the number of nodes used in the y direction and i, j, k are indices ranging
     * from 1 to nt where nt is the total number of nodes in the 3D grid, i.e.: nt = nx*ny*nz
     * with nz the number of nodes in the z direction.
     *
     * input     domain_size
     * input     grid_size
     * input     boundary_conditions
     * input     time_dep_input
     * input     gammas
     * input     source()
     * output    grid_coordinates
     * output    T
     */

    double deltax, deltay, deltaz;
    double epsilon, delold, delnew, pAp, error;
    double alpha, B;
    double **A, **L, *Ap, *y, *z, *p, *x, *xo, *r;
    int nn, nt, i, j, k, it, imax;
    double dt;
    int tt;

    /*Timing solver*/
    clock_t begin, end;
    double time_spent;
    begin = clock();

    /*Initializing parameters*/
    imax     = 5000;         //Maximum iterations ICCG
    error     = 1e-30;        //Tolerance

    nt = grid_size.nx*grid_size.ny*grid_size.nz;
    dt = (double) (time_dep_input.tf - time_dep_input.ti)/time_dep_input.timesteps;

    /*Allocating memory for computation*/
    A         = matrix2D(nt+1,4+1);
    L         = matrix2D(nt+1,4+1);
    y         = matrix1D(nt+1);
    z         = matrix1D(nt+1);
    p         = matrix1D(nt+1);
    Ap        = matrix1D(nt+1);
    x         = matrix1D(nt+1);
    xo        = matrix1D(nt+1);
    r         = matrix1D(nt+1);

    /*Initializing temperature fields*/
    for (j = 1; j <= nt; j++)
    {
        x[j]  = time_dep_input.Tinitial;
    }

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

    time_dep_input.t = time_dep_input.ti;
    tt = 0;
    do
    {
        /*Generate coefficient matrix*/
        generate_coefficient_matrix(domain_size,
                                    grid_size,
                                    boundary_conditions,
                                    time_dep_input,
                                    gammas,
                                    source,
                                    x,
                                    grid_coordinates,
                                    r,
                                    A);

        /* ICCG preconditioning */
        // Incomplete Cholesky factorization of coefficient matrix A
        incomplete_cholesky_factorization(grid_size, A, L);

        //Solving Ly=r (y = L'z)
        Ly_solver(grid_size, L, r, y);

        //Solving L'z=y
        LTz_solver(grid_size, L, y, z);

        /* Solver iterations */
        //epsilon = r'*r
        dot_product(r, r, nt, &epsilon);

        //Setting p = z
        for (j=1;j<=nt;j++)
            p[j]=z[j];

        it = 0;
        do
        {
            //Calculating Ap = A*p
            calculate_Ap(grid_size, A, p, Ap);

            //delold = r'*z
            dot_product(r, z, nt, &delold);

            //pAp = p'*Ap
            dot_product(p, Ap, nt, &pAp);

            alpha = delold/pAp;

            //x = x+alpha*p
            vector_addition(x, 1.0,    p, alpha, nt, x);

            //r = r-alpha*Ap
            vector_addition(r, 1.0,    Ap, -alpha, nt, r);

            //Solving Ly=r (y = L'z)
            Ly_solver(grid_size,L,r,y);

            //Solving L'z=y
            LTz_solver(grid_size,L,y,z);

            //delnew = r'*z
            dot_product(r, z, nt, &delnew);

            B = delnew/delold;

            //p = z + B*p;
            vector_addition(z, 1.0,    p, B, nt, p);

            //epsilon = r'*r
            dot_product(r, r, nt, &epsilon);

            //calculating error
            epsilon = sqrt(epsilon/nt);
            it = it + 1;

        }while (it < imax && epsilon > error);

        //Setting xo = x
        for (j=1;j<=nt;j++)
            xo[j]=x[j];

        time_dep_input.t = time_dep_input.t + dt;
        tt++;

    }while (tt < time_dep_input.timesteps);

    /* Processing results */
    for (i = 1; i <= grid_size.nx ; i++)
    for (j = 1; j <= grid_size.ny ; j++)
    for (k = 1; k <= grid_size.nz; k++)
    {
        nn = i + (j-1)*grid_size.nx + (k-1)*grid_size.nx*grid_size.ny;
        T[i][j][k] = x[nn];
    }

    /* Freeing memory */
    free_memory_1D(Ap);
    free_memory_1D(y);
    free_memory_1D(z);
    free_memory_1D(p);
    free_memory_1D(x);
    free_memory_1D(xo);
    free_memory_1D(r);
    free_memory_2D(A, nt);
    free_memory_2D(L, nt);

    /* Print out some results */
    end = clock();
    time_spent = (double) (end - begin)/CLOCKS_PER_SEC;

    printf("error: %E\n", epsilon);
    printf("iterations: %d\n", it);
    printf("running time: %f\n", time_spent);

}
