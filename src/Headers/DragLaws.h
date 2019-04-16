#include <aligned_new>
//=================================================================================================
//  DragLaws.h
//  Contains the implementation for the dust drag law classes.
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
//=================================================================================================
#ifndef __DRAG_LAWS_H__
#define __DRAG_LAWS_H__



//=================================================================================================
//  Class FixedDrag
/// \brief   FixedDrag class definition.
/// \details Constant Stopping time t_s = 1/K_d
/// \author  R. A. Booth
/// \date    17/10/2015
//=================================================================================================
class FixedDrag
{
public:
	FixedDrag(const FLOAT DragCoeff)
	: _t_stop(1/DragCoeff)
	{ } ;

	FLOAT operator()(FLOAT, FLOAT, FLOAT) const
	{ return _t_stop ; }
private:
	FLOAT _t_stop ;
};

//=================================================================================================
//  Class DensityDrag
/// \brief   DensityDrag class definition.
/// \details Stopping time t_s = 1 / (K_d * rho)
/// \author  R. A. Booth
/// \date    17/10/2015
//=================================================================================================
class DensityDrag
{
public:
	DensityDrag(const FLOAT DragCoeff)
	: _K_d(DragCoeff)
    { } ;

	FLOAT operator()(FLOAT grho, FLOAT drho, FLOAT) const
	{ return 1 / ((grho + drho) * _K_d) ; }
private:
	FLOAT _K_d ;
};

//=================================================================================================
//  Class EpsteinDrag
/// \brief   EpsteinDrag class definition.
/// \details Stopping time t_s = 1 / (K_d * rho * cs)
/// \author  R. A. Booth
/// \date    17/10/2015
//=================================================================================================
// EpsteinNorm = 3 sqrt(pi/8) / 4
static const double EpsteinNorm = 0.4699928014933126;
class EpsteinDrag
{
public:
	EpsteinDrag(const FLOAT DragCoeff)
	: _K_d(DragCoeff)
	{ } ;
	FLOAT operator()(FLOAT grho, FLOAT drho, FLOAT gsound) const
		{ return EpsteinNorm / ((grho + drho) * gsound * _K_d) ; }
private:
	FLOAT _K_d ;
};


//=================================================================================================
//  Class LP12_Drag
/// \brief   LP12_Drag class definition.
/// \details Stopping time using the drag coefficient from Laibe & Price (2012),
///          t_s = rho_g rho_d / K_D (rho_d + rho_g))
/// \author  R. A. Booth
/// \date    9/11/2015
//=================================================================================================
class LP12_Drag
{
public:
	LP12_Drag(const FLOAT DragCoeff)
	: _K_d(DragCoeff)
	{ } ;
	FLOAT operator()(FLOAT grho, FLOAT drho, FLOAT) const
		{ return drho * grho / ((grho + drho)  * _K_d) ; }
private:
	FLOAT _K_d ;
};


#endif//__DRAG_LAWS_H__
