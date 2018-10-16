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

typedef int bool;

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

typedef struct time_dep_input_t {
    int timesteps;
    double ti;
    double tf;
    double rho;
    double Cp;
    double Tinitial;
    double t;
}time_dep_input_t;

typedef struct gammas_t {
    double gammax;
    double gammay;
    double gammaz;
} gammas_t;

typedef enum boundary_type_t {
    NEUMANN = 0,
    DIRICHLET = 1
} boundary_type_t;

typedef struct boundary_type_faces_t {
    boundary_type_t west_boundary;
    boundary_type_t east_boundary;
    boundary_type_t south_boundary;
    boundary_type_t north_boundary;
    boundary_type_t bottom_boundary;
    boundary_type_t top_boundary;
} boundary_type_faces_t;

typedef struct fixed_boundary_funcs_t {
    double (*fixed_boundary_west) (double y, double z, double t);
    double (*fixed_boundary_east) (double y, double z, double t);
    double (*fixed_boundary_south) (double x,double z, double t);
    double (*fixed_boundary_north) (double x,double z, double t);
    double (*fixed_boundary_bottom) (double x,double y, double t);
    double (*fixed_boundary_top) (double x,double y, double t);
} fixed_boundary_funcs_t;

typedef struct flux_boundary_funcs_t {
    double (*flux_boundary_west) (double y, double z, double t);
    double (*flux_boundary_east) (double y, double z, double t);
    double (*flux_boundary_south) (double x,double z, double t);
    double (*flux_boundary_north) (double x,double z, double t);
    double (*flux_boundary_bottom) (double x,double y, double t);
    double (*flux_boundary_top) (double x,double y, double t);
} flux_boundary_funcs_t;

typedef struct boundary_conditions_t {
    boundary_type_faces_t boundary_type_faces;
    fixed_boundary_funcs_t fixed_boundary_funcs;
    flux_boundary_funcs_t flux_boundary_funcs;
} boundary_conditions_t;
#endif /* USER_TYPES_H_ */
