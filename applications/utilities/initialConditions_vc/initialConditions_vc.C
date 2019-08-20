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
    initialConditions_vc

Description
    Generates non-standard initial conditions for test cases.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Point mesh
    pointMesh pMesh(mesh);

    // Read mechanical properties dictionary
    IOdictionary mechanicalProperties
    (
        IOobject
        (
            "mechanicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Read run parameters dictionary
    IOdictionary runParameters
    (
        IOobject
        (
            "runParameters",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Test case name
    const word& tutorial(runParameters.lookup("tutorial"));

    // Read density
    const dimensionedScalar& rho(mechanicalProperties.lookup("rho"));

    // Nodal coordinates
    const vectorField& X = mesh.points();

    // Read linear momentum field
    pointVectorField lm
    (
        IOobject
        (
            "lm",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh,
        dimensionedVector("lm", dimensionSet(1,-2,-1,0,0,0,0), vector::zero)
    );

    // Non uniform angular velocity initialsed testcases
    if (tutorial == "twistingColumn")
    {
        pointVectorField omega
        (
            IOobject("omega", mesh),
            pMesh,
            dimensionedVector(runParameters.lookup("initialAngularVelocity"))
        );

        const scalar& PI = Foam::constant::mathematical::pi;
        dimensionedScalar height("height", dimensionSet(0,1,0,0,0,0,0), 6.0);

        forAll(mesh.points(), nodeID)
        {
            lm[nodeID] =
                rho.value()*(omega[nodeID]
               *Foam::sin(PI*X[nodeID].y()/(2*height.value()))) ^ X[nodeID];
        }
    }

    lm.write();

    Info<< "\n end\n";

    return 0;
}