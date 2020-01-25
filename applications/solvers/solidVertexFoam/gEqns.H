// Update coordinate
xe = inter.pointToEdge(x);
x += deltaT*(lm/rho);

// Compute right hand sides
rhsLm.primitiveFieldRef() = op.surfaceIntegrate(tC*magSf);
rhsF.primitiveFieldRef() = op.surfaceIntegrate((lmC/rho.value())*Sf);

forAll(mesh.boundary(), patch)
{
    forAll(mesh.boundaryMesh()[patch], facei)
    {
        const label& face = mesh.boundaryMesh()[patch].start() + facei;

        forAll(mesh.faces()[face], nodei)
        {
            const label& node = mesh.faces()[face][nodei];
            label nodeB = -1;
            label nodeC = -1;

            if (nodei == 0)
            {
                nodeB = mesh.faces()[face][1];
                nodeC = mesh.faces()[face][2];
            }
            else if (nodei == 1)
            {
                nodeB = mesh.faces()[face][2];
                nodeC = mesh.faces()[face][0];
            }
            else if (nodei == 2)
            {
                nodeB = mesh.faces()[face][0];
                nodeC = mesh.faces()[face][1];
            }

            vector t1 = P[node] & N_p.boundaryField()[patch][facei];
            vector t2 = P[nodeB] & N_p.boundaryField()[patch][facei];
            vector t3 = P[nodeC] & N_p.boundaryField()[patch][facei];

            vector lm1 = lm[node];
            vector lm2 = lm[nodeB];
            vector lm3 = lm[nodeC];

            if (patch == freePatchID)
            {
                t1 = vector::zero;
                t2 = vector::zero;
                t3 = vector::zero;
            }
            if (patch == tractionPatchID)
            {
                t1 = traction;
                t2 = traction;
                t3 = traction;
            }

            vector tC = 0.125*(6.0*t1 + t2 + t3);
            vector lmC = 0.125*(6.0*lm1 + lm2 + lm3);

            rhsLm[node] += tC*(magSf_p.boundaryField()[patch][facei]/3.0);
            rhsF[node] +=
                (lmC/rho.value())*(Sf_p.boundaryField()[patch][facei]/3.0);
        }
    }
}

forAll(mesh.points(), node)
{
    rhsLm[node] = rhsLm[node]/V[node];
    rhsF[node] = rhsF[node]/V[node];
}

if (angularMomentumConservation == "yes")
{
    rhsAm.primitiveFieldRef() = op.surfaceIntegrate((xe^tC)*magSf);
    forAll(mesh.points(), node)
    {
        rhsAm[node] = rhsAm[node]/V[node];
    }
    am.AMconservation(rhsLm, rhsLm1, rhsAm, xe, i, V);
}

lm += deltaT*rhsLm;
F += deltaT*rhsF;

lm.correctBoundaryConditions();