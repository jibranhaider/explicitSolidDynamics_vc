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

#include "gradientSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gradientSchemes, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

gradientSchemes::gradientSchemes
(
	const fvMesh& vm,
    const vectorList& Sf,
    const scalarList& V
)
:
	mesh_(vm),
    pMesh_(mesh_),
    edges_(mesh_.edges()),
	X_(mesh_.points()),
    Xown_(edges_.size()),
    Xnei_(edges_.size()),

    V_(V),
    Sf_(Sf),
    Xe_(edges_.size())
{
    // Compute edge center coordinates
    forAll(edges_, edge)
    {
        const label& own = edges_[edge][0];
        const label& nei = edges_[edge][1];

        Xown_[edge] = mesh_.points()[own];
        Xnei_[edge] = mesh_.points()[nei];
        Xe_[edge] = (mesh_.points()[own] + mesh_.points()[nei])/2.0;
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

gradientSchemes::~gradientSchemes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vectorField gradientSchemes::gradient
(
    const GeometricField<scalar, pointPatchField, pointMesh>& U
)   const
{
    vectorField Ugrad(mesh_.points().size(), vector::zero);

    forAll(edges_, edge)
    {
        label own = edges_[edge][0];
        label nei = edges_[edge][1];

        const vector& grad = 0.5*(U[own] + U[nei])*Sf_[edge];
        Ugrad[own] += grad;
        Ugrad[nei] -= grad;
    }

    forAll(mesh_.boundary(), patch)
    {
        forAll(mesh_.boundaryMesh()[patch], facei)
        {
            const label& face = mesh_.boundaryMesh()[patch].start() + facei;

            forAll(mesh_.faces()[face], nodei)
            {
                const label& node = mesh_.faces()[face][nodei];
                label nodeB = -1;
                label nodeC = -1;

                if (nodei == 0)
                {
                    nodeB = mesh_.faces()[face][1];
                    nodeC = mesh_.faces()[face][2];
                }
                else if(nodei == 1)
                {
                    nodeB = mesh_.faces()[face][2];
                    nodeC = mesh_.faces()[face][0];
                }
                else if(nodei == 2)
                {
                    nodeB = mesh_.faces()[face][0];
                    nodeC = mesh_.faces()[face][1];
                }

                const scalar& lmC = 6*U[node] + U[nodeB] + U[nodeC];
                const vector grad =
                    (1.0/24.0) * (lmC*mesh_.Sf().boundaryField()[patch][facei]);
                Ugrad[node] += grad;
            }
        }
    }

    return (Ugrad/V_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
tensorField gradientSchemes::gradient
(
   	const GeometricField<vector, pointPatchField, pointMesh>& U
) 	const
{
    tensorField Ugrad(mesh_.points().size(), tensor::zero);

    forAll(edges_, edge)
    {
        label own = edges_[edge][0];
        label nei = edges_[edge][1];

        const tensor& grad = 0.5*(U[own] + U[nei])*Sf_[edge];
        Ugrad[own] += grad;
        Ugrad[nei] -= grad;
    }

    forAll(mesh_.boundary(), patch)
    {
        forAll(mesh_.boundaryMesh()[patch], facei)
        {
            const label& face = mesh_.boundaryMesh()[patch].start() + facei;

            forAll(mesh_.faces()[face], nodei)
            {
                const label& node = mesh_.faces()[face][nodei];
                label nodeB = -1;
                label nodeC = -1;

                if (nodei == 0)
                {
                    nodeB = mesh_.faces()[face][1];
                    nodeC = mesh_.faces()[face][2];
                }
                else if(nodei == 1)
                {
                    nodeB = mesh_.faces()[face][2];
                    nodeC = mesh_.faces()[face][0];
                }
                else if(nodei == 2)
                {
                    nodeB = mesh_.faces()[face][0];
                    nodeC = mesh_.faces()[face][1];
                }

                const vector& lmC = 6*U[node] + U[nodeB] + U[nodeC];
                const tensor grad =
                    (1.0/24.0) * (lmC*mesh_.Sf().boundaryField()[patch][facei]);
                Ugrad[node] += grad;
            }
        }
    }

    return (Ugrad/V_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void gradientSchemes::reconstruct
(
    GeometricField<scalar, pointPatchField, pointMesh>& U,
    vectorField& Ugrad,
    scalarList& Uown,
    scalarList& Unei
)
{
    forAll(edges_, edge)
    {
        const label& own = edges_[edge][0];
        const label& nei = edges_[edge][1];

        Uown[edge] = U[own] + (Ugrad[own] & (Xe_[edge] - Xown_[edge]));
        Unei[edge] = U[nei] + (Ugrad[nei] & (Xe_[edge] - Xnei_[edge]));
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    GeometricField<vector, pointPatchField, pointMesh>& U,
    tensorField& Ugrad,
    vectorList& Uown,
    vectorList& Unei
)
{
    forAll(edges_, edge)
    {
        const label& own = edges_[edge][0];
        const label& nei = edges_[edge][1];

        Uown[edge] = U[own] + (Ugrad[own] & (Xe_[edge] - Xown_[edge]));
        Unei[edge] = U[nei] + (Ugrad[nei] & (Xe_[edge] - Xnei_[edge]));
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //