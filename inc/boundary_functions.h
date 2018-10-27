/*
 * boundary_functions.h
 *
 *  Created on: Oct 25, 2018
 *      Author: derek
 */

#ifndef BOUNDARY_FUNCTIONS_H_
#define BOUNDARY_FUNCTIONS_H_


double boundary_west(double y, double z, double t);
double boundary_east(double y, double z, double t);
double boundary_south(double x,double z, double t);
double boundary_north(double x,double z, double t);
double boundary_bottom(double x,double y, double t);
double boundary_top(double x,double y, double t);
double source_equation(double x, double y, double z, double t);

#endif /* BOUNDARY_FUNCTIONS_H_ */
