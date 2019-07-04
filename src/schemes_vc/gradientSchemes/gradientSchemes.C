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
    // Edge centers
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

// pointTensorField gradientSchemes::gradient
tensorField gradientSchemes::gradient
(
   	const GeometricField<vector, pointPatchField, pointMesh>& U
) 	const
{

    // tmp<GeometricField<tensor, pointPatchField, pointMesh> > ttf
    // (
    //     new GeometricField<tensor, pointPatchField, pointMesh>
    //     (
    //         IOobject
    //         (
    //             "gradient("+U.name()+')',
    //             mesh_
    //         ),
    //         pMesh_,
    //         dimensioned<tensor>
    //         (
    //             "0",
    //             U.dimensions()/dimLength,
    //             pTraits<tensor>::zero
    //         )
    //     )
    // );

    // GeometricField<tensor, pointPatchField, pointMesh> Ugrad = ttf();


    tensorField Ugrad(mesh_.points().size(), tensor::zero);

    forAll(edges_, edgeID)
    {
        label ownID = edges_[edgeID][0];
        label neiID = edges_[edgeID][1];

        const tensor& grad = 0.5 * (U[ownID]+U[neiID]) * Sf_[edgeID];
        Ugrad[ownID] += grad;
        Ugrad[neiID] -= grad;
    }

    forAll(mesh_.boundary(), patchID)
    {
        forAll(mesh_.boundaryMesh()[patchID], face)
        {
            const label& faceID = mesh_.boundaryMesh()[patchID].start() + face;

            forAll(mesh_.faces()[faceID], node)
            {
                const label& nodeID = mesh_.faces()[faceID][node];
                label nodeB = -1;
                label nodeC = -1;

                if (node == 0)
                {
                    nodeB = mesh_.faces()[faceID][1];
                    nodeC = mesh_.faces()[faceID][2];
                }
                else if(node == 1)
                {
                    nodeB = mesh_.faces()[faceID][2];
                    nodeC = mesh_.faces()[faceID][0];
                }
                else if(node == 2)
                {
                    nodeB = mesh_.faces()[faceID][0];
                    nodeC = mesh_.faces()[faceID][1];
                }

                const vector& lmC = 6*U[nodeID] + U[nodeB] + U[nodeC];
                const tensor grad = (1.0/24.0) * (lmC * mesh_.Sf().boundaryField()[patchID][face]);
                Ugrad[nodeID] += grad;
            }
        }
    }

    // ttf.clear();

    return (Ugrad/V_);





 //    tmp<GeometricField<vector, pointPatchField, pointMesh> > tvf
 //    (
 //        new GeometricField<vector, pointPatchField, pointMesh>
 //        (
 //            IOobject
 //            (
 //                "dummy",
 //                mesh_
 //            ),
 //            U.mesh(),
 //            dimensioned<vector>("dummy", U.dimensions()/dimLength, pTraits<vector>::zero)
 //        )
 //    );

	// GeometricField<vector, pointPatchField, pointMesh> UgradX = tvf();
	// GeometricField<vector, pointPatchField, pointMesh> UgradY = tvf();
	// GeometricField<vector, pointPatchField, pointMesh> UgradZ = tvf();

	// UgradX = gradientSchemes::gradient(U.component(0));
	// UgradY = gradientSchemes::gradient(U.component(1));
	// UgradZ = gradientSchemes::gradient(U.component(2));

 //    tmp<GeometricField<tensor, pointPatchField, pointMesh> > ttf
 //    (
 //        new GeometricField<tensor, pointPatchField, pointMesh>
 //        (
 //            IOobject
 //            (
 //                "dummy",
 //                mesh_
 //            ),
 //            U.mesh(),
 //            dimensioned<tensor>("dummy", U.dimensions()/dimLength, pTraits<tensor>::zero)
 //        )
 //    );

	// GeometricField<tensor, pointPatchField, pointMesh> Ugrad = ttf();


	// forAll(mesh_.points(), nodeID)
	// {
	// 	Ugrad[nodeID] = tensor(UgradX[nodeID], UgradY[nodeID], UgradZ[nodeID]);
	// }

	// if (Pstream::parRun())
	// {
	// 	Ugrad.correctBoundaryConditions();
	// }

	// tvf.clear();
	// ttf.clear();

	// return Ugrad;

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