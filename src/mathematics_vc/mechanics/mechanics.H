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

Class
    Foam::mechanics

Description
    Evaluation of quantities widely used in continuum mechanics.

SourceFiles
    mechanics.C

\*---------------------------------------------------------------------------*/

#ifndef mechanics_H
#define mechanics_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "fvc.H"
#include "operations.H"
#include "interpolationSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

class mechanics
{
    // Private data

        // Mesh
        const fvMesh& mesh_;

        // Point mesh
        pointMesh pMesh_;

        // Eigen-values
        vector  eigVal_;

        // Eigen-vectors
        tensor  eigVec_;

        // Operations class
        operations op;

        // Interpolation class
        interpolationSchemes inter;

        // Dual material area vector
        const vectorList& Sf_;

        // Dual material normal
        vectorList N_;

private:

    //- Disallow default bitwise copy construct
        mechanics(const mechanics&);

    //- Disallow default bitwise assignment
        void operator=(const mechanics&);

public:

    // Declare name of the class and its debug switch
        TypeName("mechanics");

    // Constructors
        mechanics
        (
            const fvMesh&,
            const vectorList&
        );


    // Destructor
        virtual ~mechanics();

    // Member functions

        // Computes spatial surface normals for the dual mesh
        // vectorList spatialNormal();
        vectorField spatialNormal();

        // Computes stretch
        pointScalarField stretch();
        // scalarField stretch();

        // Compute stabilisation
        tensorField Smatrix
        (
            const GeometricField<scalar, pointPatchField, pointMesh>&,
            const GeometricField<scalar, pointPatchField, pointMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
