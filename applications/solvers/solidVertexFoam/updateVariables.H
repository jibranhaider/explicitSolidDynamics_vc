// Strain measures
H = det(F)*op.invT(F);
J = det(F);
stretch = mech.stretch();

// Wave speeds
forAll(mesh.points(), node)
{
    Up_time[node] = model.Up().value()/stretch[node];
}

// Constitutive model
model.correct();
P = model.piola();

// Linear reconstruction
lmGrad = grad.gradient(lm);
grad.reconstruct(lm, lmGrad, lmM, lmP);

// Contact fluxes
lmC = inter.pointToEdge(lm);
tC = (inter.pointToEdge(P) & N) + (0.5*mech.Smatrix(Up, Us) & (lmP-lmM));