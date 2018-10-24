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

    clock_t begin, end;
    double time_spent;
    kershaw_algorithm_data_t* kershaw_data = NULL;


    /*Allocating memory*/
    kershaw_data = allocate_kershaw_data(grid_size);


    /*Initialize data*/
    initialize_temperature_field(grid_size,
                                 time_dep_input,
                                 kershaw_data->x);

    generate_grid_coordinates(domain_size,
                              grid_size,
                              grid_coordinates);

    initialize_time_data(time_dep_input);


    /*Execute kershaw algorithm*/
    begin = clock();

    do
    {
        generate_coefficient_matrix(domain_size,
                                    grid_size,
                                    boundary_conditions,
                                    time_dep_input,
                                    gammas,
                                    source,
                                    grid_coordinates,
                                    kershaw_data);

        preconditioning(grid_size,
                        kershaw_data);

        execute_kershaw_algorithm(grid_size,
                                  kershaw_data);

        time_dep_input.t = time_dep_input.t + time_dep_input.dt;
        time_dep_input.current_timestep++;

    }while (time_dep_input.current_timestep < time_dep_input.timesteps);


    /*process results*/
    processing_results(grid_size,
                       kershaw_data,
                       T);


    /*deallocate data*/
    free_kershaw_data(kershaw_data, grid_size);


    /*print some results*/
    end = clock();
    time_spent = (end - begin)/CLOCKS_PER_SEC;

    printf("error: %E\n", kershaw_data->epsilon);
    printf("iterations: %d\n", kershaw_data->iterations);
    printf("running time: %f\n", time_spent);

}
