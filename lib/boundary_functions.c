/*
 * boundary_functions.c
 *
 *  Created on: Oct 25, 2018
 *      Author: derek
 */


/*-----------------------------------------------------------------------------------------------*/
double boundary_west(double y, double z, double t)
{
    /*
     * Specify the west face temperature or flux distribution for the transient heat conduction
     * equation.
     * If the boundary is of type DIRICHLET, specify the temperature distribution. If the boundary
     * is of type NEUMANN specify the flux distribution.
     *
     * The distribution or flux is of the form f(y,z,t).
     *
     * input    y
     * input    z
     * input    t
     *
     * return   temperature_west
     */

    double temperature_west = 300.0;

    return temperature_west;

}


/*-----------------------------------------------------------------------------------------------*/
double boundary_east(double y, double z, double t)
{
    /*
     * Specify the east face temperature  or flux distribution for the transient heat conduction
     * equation.
     * If the boundary is of type DIRICHLET, specify the temperature distribution. If the boundary
     * is of type NEUMANN specify the flux distribution.
     *
     * The distribution or flux is of the form f(y,z,t).
     *
     * input    y
     * input    z
     * input    t
     *
     * return   temperature_east
     */

    double temperature_east = 800.0;

    return temperature_east;

}


/*-----------------------------------------------------------------------------------------------*/
double boundary_south(double x,double z, double t)
{
    /*
     * Specify the south face temperature or flux distribution for the transient heat conduction
     * equation.
     * If the boundary is of type DIRICHLET, specify the temperature distribution. If the boundary
     * is of type NEUMANN specify the flux distribution.
     *
     * The distribution or flux is of the form f(x,z,t).
     *
     * input    x
     * input    z
     * input    t
     *
     * return   temperature_south
     */

    double temperature_south = 300.0;

    return temperature_south;

}


/*-----------------------------------------------------------------------------------------------*/
double boundary_north(double x,double z, double t)
{
    /*
     * Specify the north face temperature or flux distribution for the transient heat conduction
     * equation.
     * If the boundary is of type DIRICHLET, specify the temperature distribution. If the boundary
     * is of type NEUMANN specify the flux distribution.
     *
     * The distribution or flux is of the form f(x,z,t).
     *
     * input    x
     * input    z
     * input    t
     *
     * return   temperature_north
     */

    double temperature_north = 800.0;

    return temperature_north;

}


/*-----------------------------------------------------------------------------------------------*/
double boundary_bottom(double x,double y, double t)
{
    /*
     * Specify the bottom face temperature or flux distribution for the transient heat conduction
     * equation.
     * If the boundary is of type DIRICHLET, specify the temperature distribution. If the boundary
     * is of type NEUMANN specify the flux distribution.
     *
     * The distribution or flux is of the form f(x,y,t).
     *
     * input    x
     * input    y
     * input    t
     *
     * return   temperature_bottom
     */

    double temperature_bottom = 800.0;

    return temperature_bottom;

}


/*-----------------------------------------------------------------------------------------------*/
double boundary_top(double x,double y, double t)
{
    /*
     * Specify the top face temperature or flux distribution for the transient heat conduction
     * equation.
     * If the boundary is of type DIRICHLET, specify the temperature distribution. If the boundary
     * is of type NEUMANN specify the flux distribution.
     *
     * The distribution or flux is of the form f(x,y,t).
     *
     * input    x
     * input    y
     * input    t
     *
     * return   temperature_top
     */

    double temperature_top = 300.0;

    return temperature_top;

}


/*-----------------------------------------------------------------------------------------------*/
double source_equation(double x, double y, double z, double t)
{
    /*
     * Specify the source equation for the transient heat conduction equation:
     *
     * gammax*d2T/dx2 + gammay*d2T/dy2 + gammaz*d2T/dz2 + q(x,y,z,t) = rho*Cp*dT/dt
     *
     * The source equation is of the form q(x,y,z,t).
     *
     * input    x
     * input    y
     * input    z
     * input    t
     *
     * return   q
     */

    double q = 0.0;

    return q;

}
