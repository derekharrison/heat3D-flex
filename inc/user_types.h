/*
 * user_types.h
 *
 *  Created on: Oct 11, 2018
 *      Author: derek
 */

#ifndef USER_TYPES_H_
#define USER_TYPES_H_


#define TRUE 1
#define FALSE 0
#define N_BOUNDARIES 6

typedef int bool;
typedef double(*func_pointer)(double, double, double);

typedef struct domain_size_t {
    double Lx;
    double Ly;
    double Lz;
} domain_size_t;

typedef struct fixed_boundaries_t {
    double Tw;
    double Te;
    double Ts;
    double Tn;
    double Tb;
    double Tt;
} fixed_boundaries_t;

typedef struct grid_size_t {
    int nx;
    int ny;
    int nz;
} grid_size_t;

typedef struct grid_coordinates_t {
    double*** X;
    double*** Y;
    double*** Z;
} grid_coordinates_t;

typedef struct time_data_t {
    int timesteps;
    int current_timestep;
    double ti;
    double tf;
    double Tinitial;
    double t;
    double dt;
} time_data_t;

typedef struct conductivity_t {
    double gammax;
    double gammay;
    double gammaz;
} conductivity_t;

typedef struct physical_paramaters_t {
    conductivity_t conductivity;
    double rho;
    double Cp;
} physical_paramaters_t;

typedef enum boundary_type_t {
    NEUMANN = 0,
    DIRICHLET = 1
} boundary_type_t;

typedef enum boundary_t {
    WEST = 0,
    EAST = 1,
    SOUTH = 2,
    NORTH = 3,
    BOTTOM = 4,
    TOP = 5
} boundary_t;

typedef struct boundary_type_faces_t {
    boundary_type_t west_boundary;
    boundary_type_t east_boundary;
    boundary_type_t south_boundary;
    boundary_type_t north_boundary;
    boundary_type_t bottom_boundary;
    boundary_type_t top_boundary;
} boundary_type_faces_t;

typedef struct boundary_funcs_t {
    double (*boundary_west) (double y, double z, double t);
    double (*boundary_east) (double y, double z, double t);
    double (*boundary_south) (double x,double z, double t);
    double (*boundary_north) (double x,double z, double t);
    double (*boundary_bottom) (double x,double y, double t);
    double (*boundary_top) (double x,double y, double t);
} boundary_funcs_t;

typedef struct boundary_func_type_t {
    double (*fixed_boundary) (double x1, double x2, double t);
    double (*flux_boundary) (double x1, double x2, double t);
} boundary_func_type_t;

typedef struct boundary_conditions_t {
    boundary_type_faces_t boundary_type_faces;
    boundary_funcs_t boundary_funcs;
} boundary_conditions_t;

typedef struct kershaw_algorithm_data_t {
    double** A;
    double* r;
    double** L;
    double* y;
    double* z;
    double* Ap;
    double* p;
    double* lltr;
    double* x;
    double epsilon;
    int iterations;
} kershaw_algorithm_data_t;

#endif /* USER_TYPES_H_ */
