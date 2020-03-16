/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "angularMomentum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(angularMomentum, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

angularMomentum::angularMomentum
(
    const fvMesh& vm,
    const dictionary& dict
)
:
    mesh_(vm),
    pMesh_(mesh_),
    rho_(dict.lookup("rho"))
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

angularMomentum::~angularMomentum()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void angularMomentum::AMconservation
(
    GeometricField<vector, pointPatchField, pointMesh>& rhsLm,
    GeometricField<vector, pointPatchField, pointMesh>& rhsLm1,
    const GeometricField<vector, pointPatchField, pointMesh>& rhsAm,
    scalar RKstage,
    const GeometricField<scalar, pointPatchField, pointMesh>& V
) const
{

    const objectRegistry& db = mesh_.thisDb();
    const pointVectorField& x_ = db.lookupObject<pointVectorField>("x");
    const pointVectorField& lm_ = db.lookupObject<pointVectorField>("lm");

    const dimensionedScalar deltaT
    (
        "deltaT",
        dimensionSet(0,0,1,0,0,0,0),
        db.time().deltaTValue()
    );

    tmp<GeometricField<vector, pointPatchField, pointMesh> > tvf_x
    (
        new GeometricField<vector, pointPatchField, pointMesh>
        (
            IOobject
            (
                "xAM",
                x_.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh_,
            dimensioned<vector>("xAM", x_.dimensions(), pTraits<vector>::zero)
        )
    );

    GeometricField<vector, pointPatchField, pointMesh> xAM = tvf_x();

    tmp<GeometricField<vector, pointPatchField, pointMesh> > tvf_lm
    (
        new GeometricField<vector, pointPatchField, pointMesh>
        (
            IOobject
            (
                "lmAM",
                lm_.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh_,
            dimensioned<vector>("lmAM", lm_.dimensions(), pTraits<vector>::zero)
        )
    );

    GeometricField<vector, pointPatchField, pointMesh> lmAM = tvf_lm();

    if (RKstage == 0)
    {
        xAM.primitiveFieldRef() = x_.oldTime();
    }
    else if (RKstage == 1)
    {
        xAM.primitiveFieldRef() = x_.oldTime() + (deltaT/2.0)*(lm_.oldTime()/rho_);
        lmAM.primitiveFieldRef() = lm_.oldTime() + (deltaT*rhsLm1);
        xAM.primitiveFieldRef() = xAM + ((deltaT*(lmAM/rho_))/2.0);
    }

    tensor K_LL = tensor::zero;
    tensor K_LB = tensor::zero;
    scalar K_BB = 0.0;
    vector R_L = vector::zero;

    forAll(mesh_.points(), node)
    {
        K_LL +=
            V[node]
           *((xAM[node] & xAM[node])*tensor::I - (xAM[node]*xAM[node]));

        K_LB +=
            V[node]
           *tensor(0, -xAM[node].z(), xAM[node].y(), xAM[node].z(), 0, -xAM[node].x(), -xAM[node].y(), xAM[node].x(), 0);

        K_BB += -V[node];

        R_L +=
            (V[node]*rhsAm[node])
          + ((V[node]*rhsLm[node]) ^ xAM[node]);
    }

    if (Pstream::parRun())
    {
        reduce(K_LL, sumOp<tensor>());
        reduce(K_LB, sumOp<tensor>());
        reduce(K_BB, sumOp<scalar>());
        reduce(R_L, sumOp<vector>());
    }

    tensor LHS = K_LL - ((K_LB & K_LB)/K_BB);
    vector RHS = R_L;

    vector lambda = inv(LHS) & RHS;
    vector beta = (-K_LB & lambda)/K_BB;

    forAll(mesh_.points(), node)
    {
        rhsLm[node] = rhsLm[node] + (lambda ^ xAM[node]) + beta;
    }

    /*// Constraint to check angular momentum conservation for debugging
    vector sumAM = vector::zero;
    forAll(mesh_.points(), node){
        sumAM += (V[node]*rhsAm[node]) - V[node]*(xAM[node] ^ rhsLm[node]);
    }
    Info<< " Angular momentum constraint (RK stage " << RKstage << ") = "
        << sumAM << endl;*/

    if (RKstage == 0)
    {
        rhsLm1 = rhsLm;
    }

    tvf_x.clear();
    tvf_lm.clear();

    if (Pstream::parRun())
    {
        rhsLm.correctBoundaryConditions();
        rhsLm1.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void angularMomentum::printGlobalMomentum
(
    const GeometricField<vector, pointPatchField, pointMesh>& lm,
    const GeometricField<vector, pointPatchField, pointMesh>& x,
    const GeometricField<scalar, pointPatchField, pointMesh>& V
) const
{

    vector lmG = vector::zero;
    vector amG = vector::zero;
    scalar vol = gSum(V);

    forAll(mesh_.points(), node)
    {
        lmG += lm[node]*V[node];
        amG += V[node]*(x[node] ^ lm[node]);
    }

    if (Pstream::parRun())
    {
        reduce(lmG, sumOp<vector>());
        reduce(amG, sumOp<vector>());
    }

    Info<< " Global linear momentum = " << lmG/vol << nl
        << " Global angular momentum = " << amG/vol << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //