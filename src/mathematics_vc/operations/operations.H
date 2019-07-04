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
    Foam::operations

Description
    Definition of commonly used mathematical operations.

SourceFiles
    operations.C

\*---------------------------------------------------------------------------*/

#ifndef operations_H
#define operations_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "fvc.H"
#include "vectorList.H"
#include "tensorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;


class operations
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

private:

    //- Disallow default bitwise copy construct
        operations(const operations&);

    //- Disallow default bitwise assignment
        void operator=(const operations&);

public:

    // Declare name of the class and its debug switch
        TypeName("operations");

    // Constructors
        operations
        (
            const fvMesh&
        );


    // Destructor
        virtual ~operations();

    // Member functions

        const vector& eigenValue() const
        {
            return eigVal_;
        }

        const tensor& eigenVector() const
        {
            return eigVec_;
        }


        // Minimum edge length
        dimensionedScalar minimumEdgeLength();


        // Inverse transpose of a pointTensorField
        pointTensorField invT
        (
            const GeometricField<tensor, pointPatchField, pointMesh>& T
        ) const;

        // // Tensor cross product of a volTensorField with another volTensorField
        // volTensorField tensorProduct
        // (
        //     const GeometricField<tensor, fvPatchField, volMesh>& T1,
        //     const GeometricField<tensor, fvPatchField, volMesh>& T2
        // ) const;

        // // Tensor cross product of a tensor with another tensor
        // tensor tensorProduct
        // (
        //     const tensor& T1,
        //     const tensor& T2
        // ) const;


        // // Decompose tensor
        // void decomposeTensor
        // (
        //     const GeometricField<tensor, fvPatchField, volMesh>& T,
        //     GeometricField<vector, fvPatchField, volMesh>& Vx,
        //     GeometricField<vector, fvPatchField, volMesh>& Vy,
        //     GeometricField<vector, fvPatchField, volMesh>& Vz
        // ) const;

        // // Decompose tensor and return X row
        // pointVectorField decomposeTensorX
        // (
        //     const GeometricField<tensor, pointPatchField, pointMesh>& T
        // ) const;

        // // Decompose tensor and return Y row
        // pointVectorField decomposeTensorY
        // (
        //     const GeometricField<tensor, pointPatchField, pointMesh>& T
        // ) const;

        // // Decompose tensor and return Z row
        // pointVectorField decomposeTensorZ
        // (
        //     const GeometricField<tensor, pointPatchField, pointMesh>& T
        // ) const;

        // Computes eigen-values and eigen-vectors
        void eigenStructure(const tensor&);

        // Computes area vectors for an arbitrary face
        void areaVectors
        (
            const vectorList&, vector &, scalar&
        ) const;

        // Surface integration for a vectorList
        // vectorList surfaceIntegrate
        pointVectorField surfaceIntegrate
        (
            const vectorList&
        );

        // Surface integration for a tensorList
        // tensorList surfaceIntegrate
        pointTensorField surfaceIntegrate
        (
            const tensorList&
        ) const;
};




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
