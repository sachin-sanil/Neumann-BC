 // $Id$
 /**

 * @file ./EDIT

 * @brief Definitions of Grid_int class  
 
 * @author Mayank, Sachin, Shivang, Vinayak 
 * @author Mat. No. : 22040622 EDIT
 */
 
 
///@brief Grid_int Class Declarations in header file Grid_int.h. Grid generation and applying boundary conditions.
///	Class definitions in Grid_int.cpp


#ifndef GRID_H
#define GRID_H

#include<iostream>
#include<vector>
#include<cmath>


/// @param h= grid spacing
/// @param ngp= number of grid points in one row = (2^l+1).
/// @param level =  level of multigrid = l
/// @param frc = force vector
/// @param u_app = U approximation vector
/// @param res = residual vector




class Grid
{
double h,level,ngp;
  

public:

///@brief Constructor Declarations.
///@param = l_level Local grid level
Grid(int l_level);
std::vector<double> frc;
std::vector<double> u_app;
std::vector<double> res;

/// To access private data menbers ///
double get_hValue();
double get_ngpValue();
std::vector<double> get_u_app();

};

#endif

