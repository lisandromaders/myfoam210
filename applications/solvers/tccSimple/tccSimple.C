/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    FANZYfgm

Description
    Transient solver for laminar or turbulent flow of compressible fluids
    for HVAC and similar applications.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient simulations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicPsiThermo.H"       
#include "turbulenceModel.H"
 
#include "fanzyLookUp.H"
#include "hPsiFanzy.H"

int main(int argc, char *argv[])
{
 
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "initFGM.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
 
#   include "CourantNo.H"
 
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
 
        // Reading control parameters for the SIMPLE algorithm 
	#include "simpleControl.H"
 
        // Performing non-orthogonal correction loop
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
		#include "UEqn.H"
		#include "pvEqn.H"
		#include "ztEqn.H"
        }
 
        // Writing data
        if (runTime.write())
        { 

            #include "writeFGMfields.H"
            #include "writeThermoPropertyFields.H"

        }
    }
 
    Info<< "End\n" << endl;
 
    return 0;
}

// ************************************************************************* //
