/*
 * export_data.h
 *
 *  Created on: Sep 25, 2018
 *      Author: Derek W. Harrison
 */

#ifndef EXPORT_DATA_H_
#define EXPORT_DATA_H_

#include "user_types.h"


void export_temperature_field(char* file_name,
                              char* type_data,
                              grid_size_t grid_size,
                              double ***ptr);

void export_T_data_alongz_atxy(char* file_name,
                               grid_size_t grid_size,
                               grid_coordinates_t* grid_coordinates,
                               int node_x,
                               int node_y,
                               double ***T);

void export_T_data_alongy_atxz(char* file_name,
                               grid_size_t grid_size,
                               grid_coordinates_t* grid_coordinates,
                               int node_x,
                               int node_z,
                               double ***T);

void export_T_data_alongx_atyz(char* file_name,
                               grid_size_t grid_size,
                               grid_coordinates_t* grid_coordinates,
                               int node_y,
                               int node_z,
                               double ***T);

void write_data_to_vtk(char* file_name,
                       grid_size_t grid_size,
                       grid_coordinates_t* grid_coordinates,
                       double ***ptr);

void export_data(grid_size_t grid_size,
                 grid_coordinates_t* grid_coordinates,
                 double ***T);

#endif /* EXPORT_DATA_H_ */
