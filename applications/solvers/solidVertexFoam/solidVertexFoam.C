/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    solidVertexFoam

Description
    A large strain solid mechanics solver based on a linear momentum/strains
    mixed formualtion. An explicit Total Lagrangian formulation utilisiing
    a monolithic Total Variation Diminishing Runge-Kutta time integrator.
    A discrete angular momentum projection algorithm based on two global
    Lagrange Multipliers is added for angular momentum conservation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointFields.H"
#include "gradientSchemes.H"
#include "interpolationSchemes.H"
#include "angularMomentum.H"
#include "solidModel.H"
#include "mechanics.H"
#include "dualMesh.H"
#include "operations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readControls.H"
    #include "createFields.H"

    while (runTime.loop())
    {
        if (timeStepping == "variable")
        {
            deltaT = (cfl*h)/model.Up();
            runTime.setDeltaT(deltaT);
        }

        t += deltaT; tstep++;

        Info << "\nTime Step =" << tstep << "\ndeltaT = " << deltaT.value() << " s"
             << "\nTime = " << t.value() << " s" << endl;

        lm.oldTime();
        F.oldTime();
        x.oldTime();

        forAll(RKstage, i)
        {
            #include "gEqns.H"

            if (RKstage[i] == 0)
            {
                #include "updateVariables.H"
            }
        }

        lm = 0.5*(lm.oldTime() + lm);
        F = 0.5*(F.oldTime() + F);
        x = 0.5*(x.oldTime() + x);

        #include "updateVariables.H"

        if (runTime.outputTime())
        {
            u = x - X;
            u.write();

            p = model.pressure();
            p.write();
        }

        Info << "Percentage completed = "
             << (t.value()/runTime.endTime().value())*100 << "%" << endl;
    }

    Info << "\nExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //