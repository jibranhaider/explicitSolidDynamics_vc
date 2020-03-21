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
    Foam::gradientSchemes

Description
    Define gradient computation schemes.

SourceFiles
    gradientSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef gradientSchemes_H
#define gradientSchemes_H

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class gradientSchemes
{
    // Private data

        //- Mesh
       	const fvMesh& mesh_;

        //- Point mesh
        pointMesh pMesh_;

        //- Edge list
        const edgeList& edges_;

        //- Nodal coordinate
        const vectorField& X_;

        //- Owner node coordinate
        vectorList Xown_;

        //- Neighbour node coordinate
        vectorList Xnei_;

        //- Dual control volume
        const scalarList& V_;

        //- Dual area vector
        const vectorList& Sf_;

        //- Edge center coordinate
        vectorList Xe_;

private:

    // Private member functions

        //- Disallow default bitwise copy construct
        gradientSchemes(const gradientSchemes&);

        //- Disallow default bitwise assignment
        void operator=(const gradientSchemes&);

public:

    //- Runtime type information
  	TypeName("gradientSchemes");


    // Constructors

        //- Construct from mesh
        gradientSchemes(const fvMesh&, const vectorList&, const scalarList&);


    //- Destructor
    virtual ~gradientSchemes();


    // Member functions

        // Edit

            //- Least square gradient of a pointScalarField
            vectorField gradient
            (
                const GeometricField<scalar, pointPatchField, pointMesh>& U
            ) const;

            //- Least square gradient of a pointVectorField
            tensorField gradient
            (
                const GeometricField<vector, pointPatchField, pointMesh>& U
            ) const;

            //- Linear recontruction of a pointScalarField
            void reconstruct
            (
                GeometricField<scalar, pointPatchField, pointMesh>&,
                vectorField&,
                scalarList&,
                scalarList&
            );

            //- Linear recontruction of a pointVectorField
            void reconstruct
            (
                GeometricField<vector, pointPatchField, pointMesh>&,
                tensorField&,
                vectorList&,
                vectorList&
            );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //