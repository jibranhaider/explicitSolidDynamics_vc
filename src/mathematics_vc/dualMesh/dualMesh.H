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
    Foam::dualMesh

Description
    Evaluation of quantities widely used in continuum dualMesh.

SourceFiles
    dualMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dualMesh_H
#define dualMesh_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "fvc.H"
#include "operations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

class dualMesh
{
    // Private data

        // Mesh
        const fvMesh& mesh_;

        // Point mesh
        pointMesh pMesh_;

        //- Edge list
        const edgeList& edges_;

        // Dual control volume
        pointScalarField V_;
        // scalarList V_;

        // Edge centers
        vectorList Xe_;

        // Edge-cell-face connectivity
        List<List<labelList>> edgeCellsFaces_;

        // Interior nodes
        PackedBoolList isInteriorNode_;

        // Interior edges
        PackedBoolList isInteriorEdge_;

        // Dual material face area vectors
        vectorList Sf_;

        // Operations class
        operations op;

private:

    //- Disallow default bitwise copy construct
        dualMesh(const dualMesh&);

    //- Disallow default bitwise assignment
        void operator=(const dualMesh&);

public:

    // Declare name of the class and its debug switch
        TypeName("dualMesh");

    // Constructors
        dualMesh
        (
            const fvMesh&
        );


    // Destructor
        virtual ~dualMesh();

    // Member functions

        // Access

            const pointScalarField& volume() const
            // const scalarList& volume() const
            {
                return V_;
            }

            const vectorList& edgeCentre() const
            {
                return Xe_;
            }

            const List<List<labelList>>& connectivity() const
            {
                return edgeCellsFaces_;
            }

            const vectorList& areaVector() const
            {
                return Sf_;
            }

            const PackedBoolList& isInteriorEdge() const
            {
                return isInteriorEdge_;
            }

        // Print

            // Print centroid of primal mesh
            void printPrimalMeshCentroid() const;

            // Print centroid of dua; mesh
            void printDualMeshCentroid() const;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
