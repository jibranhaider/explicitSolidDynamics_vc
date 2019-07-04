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

#include "dualMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
defineTypeNameAndDebug(dualMesh, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //
dualMesh::dualMesh
(
    const fvMesh& vm
)
:
    mesh_(vm),
    pMesh_(mesh_),
    edges_(mesh_.edges()),

    V_
    (
        IOobject
        (
            "V_",
            mesh_
        ),
        pMesh_,
        dimensionedScalar("V_", dimVolume, 0.0)
    ),

    Xe_(edges_.size()),
    edgeCellsFaces_(edges_.size()),
    isInteriorNode_(mesh_.nPoints(), true),
    Sf_(edges_.size()),
    op(mesh_)
{

    // Dual control volume
    forAll(mesh_.cells(), cell)
    {
        const scalar& size = mesh_.cellPoints()[cell].size();
        const scalar& V_t = mesh_.V()[cell]/size;

        forAll(mesh_.cellPoints()[cell], point)
        {
            const label& node = mesh_.cellPoints()[cell][point];
            V_[node] += V_t;
        }
    }

    // Edge centers
    forAll(edges_, edge)
    {
        const label& own = edges_[edge][0];
        const label& nei = edges_[edge][1];
        Xe_[edge] = (mesh_.points()[own] + mesh_.points()[nei])/2.0;
    }

    // Edge-cell-face connectivity
    forAll(mesh_.edgeCells(), edge)
    {
        edgeCellsFaces_[edge].setSize(mesh_.edgeCells()[edge].size());

        forAll(mesh_.edgeCells()[edge], celli)
        {
            const label& cell = mesh_.edgeCells()[edge][celli];

            DynamicList<label> faces(0);

            forAll(mesh_.edgeFaces()[edge], facei)
            {
                const scalar& face = mesh_.edgeFaces()[edge][facei];

                forAll(mesh_.cells()[cell], i)
                {
                    const scalar& faceCheck = mesh_.cells()[cell][i];

                    if (face == faceCheck)
                    {
                        faces.append(face);
                    }
                }
            }

            edgeCellsFaces_[edge][celli].setSize(faces.size());

            forAll(edgeCellsFaces_[edge][celli], facei)
            {
                edgeCellsFaces_[edge][celli][facei] = faces(facei);
            }
        }
    }

    // Interior nodes
    forAll(mesh_.boundary(), patch)
    {
        forAll(mesh_.boundaryMesh()[patch].meshPoints(), nodei)
        {
            const label& node = mesh_.boundaryMesh()[patch].meshPoints()[nodei];
            isInteriorNode_.unset(node);
        }
    }

    // Area vectors
    forAll(mesh_.edgeCells(), edge)
    {
        const label& own = edges_[edge][0];
        const label& nei = edges_[edge][1];
        const vector test = mesh_.points()[nei] - mesh_.points()[own];

        forAll(mesh_.edgeCells()[edge], celli)
        {
            const label& cell = mesh_.edgeCells()[edge][celli];

            vector n = vector::zero;
            scalar a = 0.0;

            vectorList x(4);
            x[0] = mesh_.C()[cell];
            x[2] = Xe_[edge];

            const label& face1 = edgeCellsFaces_[edge][celli][0];
            const label& face2 = edgeCellsFaces_[edge][celli][1];

            if (mesh_.isInternalFace(face1))
            {
                x[1] = mesh_.Cf()[face1];
            }
            else
            {
                const label& patch = mesh_.boundaryMesh().whichPatch(face1);
                const label& facei = mesh_.boundaryMesh()[patch].whichFace(face1);
                x[1] = mesh_.Cf().boundaryField()[patch][facei];
            }

            if (mesh_.isInternalFace(face2))
            {
                x[3] = mesh_.Cf()[face2];
            }
            else
            {
                const label& patch = mesh_.boundaryMesh().whichPatch(face2);
                const label& facei = mesh_.boundaryMesh()[patch].whichFace(face2);
                x[3] = mesh_.Cf().boundaryField()[patch][facei];
            }

            op.areaVectors(x,n,a);

            if ((test & n) < 0.0)
            {
                n = -n;
            }

            Sf_[edge] += n*a;
        }
    }

}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //
dualMesh::~dualMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dualMesh::printPrimalMeshCentroid () const
{
    vector sum = vector::zero;
    scalar Vt = gSum(mesh_.V());

    forAll(mesh_.cells(), cell)
    {
        sum += mesh_.C()[cell]*mesh_.V()[cell];
    }

    if (Pstream::parRun())
    {
        reduce(sum, sumOp<vector>());
    }

    Info << "\nCentroid of primal mesh = " << sum/Vt << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void dualMesh::printDualMeshCentroid () const
{
    vector sum = vector::zero;
    scalar Vt = gSum(V_);

    forAll(mesh_.points(), node)
    {
        sum += mesh_.points()[node]*V_[node];
    }

    if (Pstream::parRun())
    {
        reduce(sum, sumOp<vector>());
    }

    Info << "\nCentroid of dual mesh = " << sum/Vt << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //