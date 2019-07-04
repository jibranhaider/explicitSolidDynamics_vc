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

\*---------------------------------------------------------------------------*/

#include "interpolationSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(interpolationSchemes, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

interpolationSchemes::interpolationSchemes(const fvMesh& vm)
:
	mesh_(vm),
    edges_(mesh_.edges())
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

interpolationSchemes::~interpolationSchemes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalarList interpolationSchemes::pointToEdge
(
    const GeometricField<scalar, pointPatchField, pointMesh>& U
) const
{

    scalarList Ue(edges_.size());

    forAll(edges_, edge)
    {
        const label& own = edges_[edge][0];
        const label& nei = edges_[edge][1];

        Ue[edge] = 0.5*(U[own] + U[nei]);
    }

    return Ue;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

vectorList interpolationSchemes::pointToEdge
(
    const GeometricField<vector, pointPatchField, pointMesh>& U
) const
{

    vectorList Ue(edges_.size());

    forAll(edges_, edge)
    {
        const label& own = edges_[edge][0];
        const label& nei = edges_[edge][1];

        Ue[edge] = 0.5*(U[own] + U[nei]);
    }

    return Ue;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tensorList interpolationSchemes::pointToEdge
(
    const GeometricField<tensor, pointPatchField, pointMesh>& U
) const
{

    tensorList Ue(edges_.size());

    forAll(edges_, edge)
    {
        const label& own = edges_[edge][0];
        const label& nei = edges_[edge][1];

        Ue[edge] = 0.5*(U[own] + U[nei]);
    }

    return Ue;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //