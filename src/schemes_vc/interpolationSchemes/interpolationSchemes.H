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
    Foam::interpolationSchemes

Description
    Define interpolation schemes.

SourceFiles
    interpolationSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef interpolationSchemes_H
#define interpolationSchemes_H

#include "volFields.H"
#include "pointFields.H"
#include "vectorList.H"
#include "tensorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class interpolationSchemes
{
    // Private data

        //- Mesh
        const fvMesh& mesh_;

        //- Edge list
        const edgeList& edges_;

private:

    // Private member functions

        //- Disallow default bitwise copy construct
        interpolationSchemes(const interpolationSchemes&);

        //- Disallow default bitwise assignment
        void operator=(const interpolationSchemes&);

public:

    //- Runtime type information
    TypeName("interpolationSchemes");

    // Constructors

        //- Construct from mesh
            explicit interpolationSchemes(const fvMesh&);

    //- Destructor
    virtual ~interpolationSchemes();


    // Member functions

        // Edit

            //- Point to edge interpolation for a scalarList
            scalarList pointToEdge
            (
                const GeometricField<scalar, pointPatchField, pointMesh>& U
            ) const;

            //- Point to edge interpolation for a vectorList
            vectorList pointToEdge
            (
                const GeometricField<vector, pointPatchField, pointMesh>& U
            ) const;

            //- Point to edge interpolation for a tensorList
            tensorList pointToEdge
            (
                const GeometricField<tensor, pointPatchField, pointMesh>& U
            ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //