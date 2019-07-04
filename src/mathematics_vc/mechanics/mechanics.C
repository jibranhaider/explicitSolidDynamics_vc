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

#include "mechanics.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
defineTypeNameAndDebug(mechanics, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //
mechanics::mechanics
(
    const fvMesh& vm,
    const vectorList& Sf
)
:
    mesh_(vm),
    pMesh_(mesh_),
    op(mesh_),
    inter(mesh_),
    Sf_(Sf),
    N_(mesh_.edges().size())
{
    N_ = Sf/mag(Sf);
}




// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //
mechanics::~mechanics()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// vectorList mechanics::spatialNormal()
vectorField mechanics::spatialNormal()
{
    const objectRegistry& db = mesh_.thisDb();
    const pointTensorField& F_ = db.lookupObject<pointTensorField>("F");

    // vectorList n(mesh_.edges().size());
    vectorField n(mesh_.edges().size());

    forAll(mesh_.edges(), edgeID)
    {
        label ownID = mesh_.edges()[edgeID][0];
        label neiID = mesh_.edges()[edgeID][1];

        const tensor& FcInv = inv((F_[ownID] + F_[neiID]) / 2.0);
        n[edgeID] = (FcInv.T() & N_[edgeID]) / (mag(FcInv.T() & N_[edgeID]));
    }

    return n;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

pointScalarField mechanics::stretch()
// scalarField mechanics::stretch()
{
    const objectRegistry& db = mesh_.thisDb();
    const pointTensorField& F_ = db.lookupObject<pointTensorField>("F");
    pointTensorField C_ = F_.T() & F_;

    tmp<GeometricField<scalar, pointPatchField, pointMesh> > tsf
    (
        new GeometricField<scalar, pointPatchField, pointMesh>
        (
            IOobject("stretch", mesh_),
            F_.mesh(),
            dimensioned<scalar>("stretch", dimless, pTraits<scalar>::one)
        )
    );

    GeometricField<scalar, pointPatchField, pointMesh> stretch = tsf();

    // tmp<scalarField> ttf
    // (
    //     new scalarField
    //     (
    //         mesh_.points().size(),
    //         0.0
    //     )
    // );

    // scalarField stretch = ttf();

    forAll(mesh_.points(), nodeID)
    {
        op.eigenStructure(C_[nodeID]);
        vector eigVal_ = op.eigenValue();
        // vector eigVal_ = Foam::eigenValues(C_[nodeID]);
        stretch[nodeID] = min(eigVal_.x(), eigVal_.y());
        stretch[nodeID] = Foam::sqrt(min( stretch[nodeID], eigVal_.z()) );
    }

    tsf.clear();

    stretch.correctBoundaryConditions();
    return stretch;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tensorField mechanics::Smatrix
(
    const GeometricField<scalar, pointPatchField, pointMesh>& Up,
    const GeometricField<scalar, pointPatchField, pointMesh>& Us
)
{

    // tmp<tensorield> S(new scalarField(mesh_.points().size(), tensor::zero));
    tensorField S(mesh_.points().size(), tensor::zero);
    vectorField n = mechanics::spatialNormal();

    S = (inter.pointToEdge(Up)*n*n) + (inter.pointToEdge(Us)*(tensor::I-(n*n)));

    return S;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //