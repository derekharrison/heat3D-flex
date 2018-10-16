/*
 * heat3D_helper.h
 *
 *  Created on: Sep 26, 2018
 *      Author: Derek W. Harrison
 */

#ifndef heat3D_HELPER_H_
#define heat3D_HELPER_H_

#include "user_types.h"


void vector_addition(double *v1,
                     double factor1,
                     double *v2,
                     double factor2,
                     int size_vector,
                     double *output_vector);

void dot_product(double *v1,
                 double *v2,
                 int size_vec,
                 double *dotproduct);

void calculate_Ap(grid_size_t grid_size,
                  double** A,
                  double* p,
                  double* Astorp);

void generate_coefficient_matrix(domain_size_t domain_size,
                                 grid_size_t grid_size,
                                 boundary_conditions_t boundary_conditions,
                                 time_dep_input_t time_dep_input,
                                 gammas_t gammas,
                                 double (*source)(double x,double y,double z,double t),
                                 double* xo,
                                 grid_coordinates_t* grid_coordinates,
                                 double* r,
                                 double** A);

void Ly_solver(grid_size_t grid_size,
               double** L,
               double* r,
               double* y);

void LTz_solver(grid_size_t grid_size,
                double** L,
                double* y,
                double* z);

void incomplete_cholesky_factorization(grid_size_t grid_size,
                                       double** A,
                                       double** L);

#endif /* heat3D_HELPER_H_ */
