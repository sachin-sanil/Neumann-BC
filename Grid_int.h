 /** @mainpage SiWiR2 Assignment 1.*/
 
 // $Id$
 /**

 * @file ./EDIT

 * @brief Definitions of Grid_int class  
 
 * @author Mayank, Sachin, Shivang, Vinayak 
 * @author Mat. No. : 22040622
 */
 
 
///@brief Grid_int Class Declarations in header file Grid_int.h. Grid generation and applying boundary conditions.
///	Class definitions in Grid_int.cpp




#ifndef GRID_INT_H
#define GRID_INT_H

#include<iostream>
#include<vector>
#include<cmath>


/// @param h= grid spacing
/// @param ngp= number of grid points in one row = (2^l+1).
/// @param level =  level of multigrid = l
/// @param X = U initial after applying BC


class Grid_int
{
private:
double h,level,ngp;  

public:

///@brief Constructor Declarations.
///@param l_level Grid level.
Grid_int(int l_level);

std::vector<double> X;


/// To access private members ///
double get_hValue();
double get_ngpValue();



void display_grid_int();

std::vector<double> get_Xvalue();

///@brief function to apply boundary conditions ///
void boundary_con(); 
std::vector<double> U_exact();

void boundary_neu();

};
#endif



