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

#include "operations.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
defineTypeNameAndDebug(operations, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //
operations::operations
(
    const fvMesh& vm
)
:
    mesh_(vm),
    pMesh_(mesh_)
{}




// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //
operations::~operations()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


dimensionedScalar operations::minimumEdgeLength()
{
    dimensionedScalar h ("h", dimensionSet(0,1,0,0,0,0,0), GREAT);

    forAll (mesh_.edges(),edgeID)
    {
        forAll (mesh_.edges()[edgeID],point)
        {
            const label& pointID = mesh_.edges()[edgeID][point];
            const label& pointID1 = mesh_.edges()[edgeID][point+1];
            const scalar& edgeLength = mag(mesh_.points()[pointID]-mesh_.points()[pointID1]);

            if (edgeLength < h.value())
            {
                h.value() = edgeLength;
            }
            break;
        }
    }

    if( Pstream::parRun() )
    {
        reduce(h.value(), minOp<scalar>());
    }

    return h;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

pointTensorField operations::invT
(
    const GeometricField<tensor, pointPatchField, pointMesh>& T
) const
{

    tmp<GeometricField<tensor, pointPatchField, pointMesh> > tsf
    (
        new GeometricField<tensor, pointPatchField, pointMesh>
        (
            IOobject("inv", mesh_),
            T.mesh(),
            dimensioned<tensor>("inv", T.dimensions(), pTraits<tensor>::one)
        )
    );

    GeometricField<tensor, pointPatchField, pointMesh> inv = tsf();

    inv = Foam::inv(T);

    tsf.clear();

    return inv.T();
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// volTensorField operations::tensorProduct
// (
//     const GeometricField<tensor, fvPatchField, volMesh>& T1,
//     const GeometricField<tensor, fvPatchField, volMesh>& T2
// ) const
// {

//     tmp<GeometricField<tensor, fvPatchField, volMesh> > tsf
//     (
//         new GeometricField<tensor, fvPatchField, volMesh>
//         (
//             IOobject("P", mesh_),
//             mesh_,
//             dimensioned<tensor>("P", T1.dimensions()*T2.dimensions(), pTraits<tensor>::one)
//         )
//     );

//     GeometricField<tensor, fvPatchField, volMesh> P = tsf();

//     forAll (mesh_.cells(), i)
//     {
//         P[i] = operations::tensorProduct(T1[i], T2[i]);
//     }

//     tsf.clear();

//     return P;
// }


// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// tensor operations::tensorProduct
// (
//     const tensor& T1, const tensor& T2
// ) const
// {
//     tensor P = tensor::zero;

//     P.xx() = (T1.yy() * T2.zz()) - (T1.yz() * T2.zy()) + (T1.zz() * T2.yy()) - (T1.zy() * T2.yz());
//     P.xy() = (T1.yz() * T2.zx()) - (T1.yx() * T2.zz()) + (T1.zx() * T2.yz()) - (T1.zz() * T2.yx());
//     P.xz() = (T1.yx() * T2.zy()) - (T1.yy() * T2.zx()) + (T1.zy() * T2.yx()) - (T1.zx() * T2.yy());

//     P.yx() = (T1.xz() * T2.zy()) - (T1.xy() * T2.zz()) + (T1.zy() * T2.xz()) - (T1.zz() * T2.xy());
//     P.yy() = (T1.zz() * T2.xx()) - (T1.zx() * T2.xz()) + (T1.xx() * T2.zz()) - (T1.xz() * T2.zx());
//     P.yz() = (T1.zx() * T2.xy()) - (T1.zy() * T2.xx()) + (T1.xy() * T2.zx()) - (T1.xx() * T2.zy());

//     P.zx() = (T1.xy() * T2.yz()) - (T1.xz() * T2.yy()) + (T1.yz() * T2.xy()) - (T1.yy() * T2.xz());
//     P.zy() = (T1.xz() * T2.yx()) - (T1.xx() * T2.yz()) + (T1.yx() * T2.xz()) - (T1.yz() * T2.xx());
//     P.zz() = (T1.xx() * T2.yy()) - (T1.xy() * T2.yx()) + (T1.yy() * T2.xx()) - (T1.yx() * T2.xy());

//     return P;
// }


// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// void operations::decomposeTensor
// (
//     const GeometricField<tensor, fvPatchField, volMesh>& T,
//     GeometricField<vector, fvPatchField, volMesh>& Vx,
//     GeometricField<vector, fvPatchField, volMesh>& Vy,
//     GeometricField<vector, fvPatchField, volMesh>& Vz
// ) const
// {
//     forAll (mesh_.cells(), cellID)
//     {
//         Vx[cellID] = vector(T[cellID].xx(), T[cellID].xy(), T[cellID].xz());
//         Vy[cellID] = vector(T[cellID].yx(), T[cellID].yy(), T[cellID].yz());
//         Vz[cellID] = vector(T[cellID].zx(), T[cellID].zy(), T[cellID].zz());
//     }

//     if( Pstream::parRun() )
//     {
//         Vx.correctBoundaryConditions();
//         Vy.correctBoundaryConditions();
//         Vz.correctBoundaryConditions();
//     }
// }


// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// pointVectorField operations::decomposeTensorX
// (
//     const GeometricField<tensor, pointPatchField, pointMesh>& T
// ) const
// {

//     tmp<GeometricField<vector, pointPatchField, pointMesh> > tvf
//     (
//         new GeometricField<vector, pointPatchField, pointMesh>
//         (
//             IOobject
//             (
//                 "dummy",
//                 mesh_
//             ),
//             T.mesh(),
//             dimensioned<vector>("dummy", T.dimensions(), pTraits<vector>::zero)
//         )
//     );

//     GeometricField<vector, pointPatchField, pointMesh> Vx = tvf();

//     forAll (mesh_.points(), nodeID)
//     {
//         Vx[nodeID] = vector(T[nodeID].xx(), T[nodeID].xy(), T[nodeID].xz());
//     }

//     tvf.clear();
//     Vx.correctBoundaryConditions();

//     return Vx;
// }


// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// pointVectorField operations::decomposeTensorY
// (
//     const GeometricField<tensor, pointPatchField, pointMesh>& T
// ) const
// {

//     tmp<GeometricField<vector, pointPatchField, pointMesh> > tvf
//     (
//         new GeometricField<vector, pointPatchField, pointMesh>
//         (
//             IOobject
//             (
//                 "dummy",
//                 mesh_
//             ),
//             T.mesh(),
//             dimensioned<vector>("dummy", T.dimensions(), pTraits<vector>::zero)
//         )
//     );

//     GeometricField<vector, pointPatchField, pointMesh> Vy = tvf();

//     forAll (mesh_.points(), nodeID)
//     {
//         Vy[nodeID] = vector(T[nodeID].yx(), T[nodeID].yy(), T[nodeID].yz());
//     }

//     tvf.clear();
//     Vy.correctBoundaryConditions();

//     return Vy;
// }


// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// pointVectorField operations::decomposeTensorZ
// (
//     const GeometricField<tensor, pointPatchField, pointMesh>& T
// ) const
// {

//     tmp<GeometricField<vector, pointPatchField, pointMesh> > tvf
//     (
//         new GeometricField<vector, pointPatchField, pointMesh>
//         (
//             IOobject
//             (
//                 "dummy",
//                 mesh_
//             ),
//             T.mesh(),
//             dimensioned<vector>("dummy", T.dimensions(), pTraits<vector>::zero)
//         )
//     );

//     GeometricField<vector, pointPatchField, pointMesh> Vz = tvf();

//     forAll (mesh_.points(), nodeID)
//     {
//         Vz[nodeID] = vector(T[nodeID].zx(), T[nodeID].zy(), T[nodeID].zz());
//     }

//     tvf.clear();
//     Vz.correctBoundaryConditions();

//     return Vz;
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void operations::eigenStructure(const tensor& ten)
{
    tensor t(ten);

    scalar it_max = 100;
    scalar it_num = 0;
    scalar rot_num = 0;
    scalar size = 3;
    scalar gapq = 0.0;

    int i,j;
    int k,l,m,p1,q = 0;

    double term, termp, termq = 0.0;
    double g,h,c,w1 = 0.0;
    double theta,thresh = 0.0;
    double s,t1,tau1 = 0.0;

    vector d( vector(t.xx(),t.yy(),t.zz()) );
    vector bw( vector(t.xx(),t.yy(),t.zz()) );
    vector zw(vector::zero);
    tensor v1(tensor::I);


    while ( it_num < it_max )
    {
        it_num += 1;

        // The convergence threshold is based on the size of the elements in the strict upper triangle of the matrix.
        thresh = 0.0;

        for ( j = 0; j < size; j++ )
        {
          for ( i = 0; i < j; i++ )
          {
            thresh = thresh + t[i+j*size] * t[i+j*size];
          }
        }

        thresh = Foam::sqrt ( thresh ) / ( 4 * size );

        if ( thresh == 0.0 )
        {
          break;
        }

        for ( p1 = 0; p1 < size; p1++ )
        {
            for ( q = p1 + 1; q < size; q++ )
            {
                gapq = 10.0 * fabs( t[p1+q*size] );
                termp = gapq + fabs( d[p1] );
                termq = gapq + fabs( d[q] );

                // Annihilate tiny offdiagonal elements.
                if ( 4 < it_num && termp == fabs( d[p1] ) && termq == fabs( d[q] ) )
                {
                  t[p1+q*size] = 0.0;
                }

                //  Otherwise, apply a rotation.
                else if ( thresh <= fabs( t[p1+q*size] ) )
                {
                    h = d[q] - d[p1];
                    term = fabs( h ) + gapq;

                    if ( term == fabs( h ) )
                    {
                        t1 = t[p1+q*size] / h;
                    }
                    else
                    {
                        theta = 0.5 * h / t[p1+q*size];
                        t1 = 1.0 / ( fabs(theta) + Foam::sqrt( 1.0 + theta * theta ) );
                        if ( theta < 0.0 )
                        {
                            t1 = - t1;
                        }
                    }

                    c = 1.0 / Foam::sqrt ( 1.0 + t1 * t1 );
                    s = t1 * c;
                    tau1 = s / ( 1.0 + c );
                    h = t1 * t[p1+q*size];

                    //  Accumulate corrections to diagonal elements.
                    zw[p1] = zw[p1] - h;
                    zw[q] = zw[q] + h;
                    d[p1] = d[p1] - h;
                    d[q] = d[q] + h;
                    t[p1+q*size] = 0.0;

                    // Rotate, using information from the upper triangle of A only.
                    for ( j = 0; j < p1; j++ )
                    {
                        g = t[j+p1*size];
                        h = t[j+q*size];
                        t[j+p1*size] = g - s * ( h + g * tau1 );
                        t[j+q*size] = h + s * ( g - h * tau1 );
                    }

                    for ( j = p1 + 1; j < q; j++ )
                    {
                        g = t[p1+j*size];
                        h = t[j+q*size];
                        t[p1+j*size] = g - s * ( h + g * tau1 );
                        t[j+q*size] = h + s * ( g - h * tau1 );
                    }

                    for ( j = q + 1; j < size; j++ )
                    {
                        g = t[p1+j*size];
                        h = t[q+j*size];
                        t[p1+j*size] = g - s * ( h + g * tau1 );
                        t[q+j*size] = h + s * ( g - h * tau1 );
                    }

                    //  Accumulate information in the eigenvector matrix.
                    for ( j = 0; j < size; j++ )
                    {
                        g = v1[j+p1*size];
                        h = v1[j+q*size];
                        v1[j+p1*size] = g - s * ( h + g * tau1 );
                        v1[j+q*size] = h + s * ( g - h * tau1 );
                    }
                    rot_num = rot_num + 1;

                }

            }
        }
    }


    // Restore upper triangle of input matrix
    for ( i = 0; i < size; i++ )
    {
      bw[i] = bw[i] + zw[i];
      d[i] = bw[i];
      zw[i] = 0.0;
    }

    //  Ascending sort the eigenvalues and eigenvectors
    for ( k = 0; k < size - 1; k++ )
    {
        m = k;
        for ( l = k + 1; l < size; l++ )
        {
            if ( d[l] < d[m] )
            {
                m = l;
            }
        }

        if ( m != k )
        {
            t1    = d[m];
            d[m] = d[k];
            d[k] = t1;
            for ( i = 0; i < size; i++ )
            {
                w1        = v1[i+m*size];
                v1[i+m*size] = v1[i+k*size];
                v1[i+k*size] = w1;
            }
        }
    }

    // Corrections for calculating inverse
    tensor sub(tensor::zero);

    for ( i=0; i<3; i++ )
    {
        if ( d[i] < SMALL )
        {
            d[i] = 1;
            sub += d[i] * vector(v1[3*i],v1[3*i+1],v1[3*i+2]) * vector(v1[3*i],v1[3*i+1],v1[3*i+2]);
        }
    }

    eigVal_ = d;
    eigVec_ = v1;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void operations::areaVectors
(
    const vectorList& x, vector& n, scalar& a
) const
{
    vectorList iso(x.size());

    vector tanVecZeta = vector::zero;
    vector tanVecEta = vector::zero;
    scalar eta_GP = 0.0;
    scalar zeta_GP = 0.0;

    if (x.size() == 4)
    {
        iso[0] = vector(-1,-1,0);
        iso[1] = vector(1,-1,0);
        iso[2] = vector(1,1,0);
        iso[3] = vector(-1,1,0);

        for (int i=0; i<4; i++)
        {
            tanVecZeta += x[i] * (iso[i].x() * (1.0+eta_GP*iso[i].y()) / 4.0);
            tanVecEta +=  x[i] * (iso[i].y() * (1.0+zeta_GP*iso[i].x()) / 4.0);
        }

        n = (tanVecZeta ^ tanVecEta)/mag(tanVecZeta ^ tanVecEta);
        a = 4.0*mag(tanVecZeta ^ tanVecEta);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// vectorList operations::surfaceIntegrate
pointVectorField operations::surfaceIntegrate
(
    const vectorList& flux
)
{

    tmp<GeometricField<vector, pointPatchField, pointMesh> > tvf
    (
        new GeometricField<vector, pointPatchField, pointMesh>
        (
            IOobject
            (
                "surfaceIntegrate()",
                mesh_
            ),
            pMesh_,
            pTraits<vector>::zero
            // dimensionedVector("test",dimensionSet(1,-2,-2,0,0,0,0),vector::zero)
            // dimensioned<vector>("0", T.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, pointPatchField, pointMesh> rhs = tvf();

    // vectorList rhs(mesh_.edges().size(), vector::zero);

    forAll(mesh_.edges(), edge)
    {
        label own = mesh_.edges()[edge][0];
        label nei = mesh_.edges()[edge][1];

        rhs[own] += flux[edge];
        rhs[nei] -= flux[edge];
    }

    tvf.clear();

    return rhs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// tensorList operations::surfaceIntegrate
pointTensorField operations::surfaceIntegrate
(
    const tensorList& flux
) const
{

    tmp<GeometricField<tensor, pointPatchField, pointMesh> > tvf
    (
        new GeometricField<tensor, pointPatchField, pointMesh>
        (
            IOobject
            (
                "surfaceIntegrate()",
                mesh_
            ),
            pMesh_,
            pTraits<tensor>::zero
            // dimensioned<tensor>("0", T.dimensions(), pTraits<tensor>::zero)
        )
    );
    GeometricField<tensor, pointPatchField, pointMesh> rhs = tvf();


    // tensorList rhs(mesh_.edges().size(), tensor::zero);

    forAll(mesh_.edges(), edge)
    {
        label own = mesh_.edges()[edge][0];
        label nei = mesh_.edges()[edge][1];

        rhs[own] += flux[edge];
        rhs[nei] -= flux[edge];
    }

    tvf.clear();

    return rhs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //