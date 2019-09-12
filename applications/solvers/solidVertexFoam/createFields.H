Info << "Creating mesh parameters ..." << endl;

// Operations class
operations op(mesh);

// Dual mesh class
dualMesh dual(mesh);

// Interpolation schemes class
interpolationSchemes inter(mesh);

// Point mesh
pointMesh pMesh(mesh);

// Dual control volume
const pointScalarField& V = dual.volume();

// Primal mesh material face area normal vector
const surfaceVectorField& Sf_p = mesh.Sf();

// Primal mesh material face area
const surfaceScalarField& magSf_p = mesh.magSf();

// Primal mesh material normal
const surfaceVectorField N_p = Sf_p/magSf_p;

// Dual area vector
const vectorList Sf = dual.areaVector();

// Dual face area
const scalarList magSf = mag(Sf);

// Dual material normals
const vectorList N = Sf/mag(Sf);

// Minimum edge length
const dimensionedScalar& h = op.minimumEdgeLength();

// Boundary patches
const polyBoundaryMesh& bm = mesh.boundaryMesh();
const label& freePatchID = bm.findPatchID("free");
const label& tractionPatchID = bm.findPatchID("traction");


Info << "Creating mesh coordinate fields ..." << endl;

// Spatial nodal coordinate
pointVectorField x
(
    IOobject("x", mesh),
    pMesh,
    dimensionedVector("x", dimensionSet(0,1,0,0,0,0,0), vector::zero)
);
x.primitiveFieldRef() = mesh.points();

// Material nodal coordinate
pointVectorField X
(
    IOobject ("X", mesh),
    x
);

// Material edge centre coordinate
const vectorList& Xe = dual.edgeCentre();

// Spatial edge centre coordiantes
vectorList xe = Xe;


Info << "Creating linear momentum fields ..." << endl;

// Nodal linear momentum
pointVectorField lm
(
    IOobject ("lm", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    pMesh
);


Info << "Creating strain measure fields ..." << endl;

// Deformation gradient tensor
pointTensorField F
(
    IOobject("F", mesh),
    pMesh,
    Foam::tensor::I
);

// RHS of F equation
pointTensorField rhsF
(
    IOobject ("rhsF", mesh),
    pMesh,
    dimensionedTensor("rhsF",F.dimensions()/dimTime,tensor::zero)
);

// Cofactor of deformation
pointTensorField H
(
    IOobject ("H", mesh),
    det(F)*op.invT(F)
);

// Jacobian of deformation
pointScalarField J
(
    IOobject("J", mesh),
    det(F)
);


Info << "Creating constitutive model ..." << endl;

// Solid model class
solidModel model(F, mechanicalProperties, mesh);
model.correct();

// Density
const dimensionedScalar& rho = model.density();

// Pressure
pointScalarField p = model.pressure();
// p.write();


// First Piola Kirchhoff stress tensor
pointTensorField P = model.piola();


Info << "Creating fields for wave speeds ..." << endl;

// Continuum mechanics class
mechanics mech(mesh, Sf);

// Longitudinal wave speed
pointScalarField Up
(
    IOobject("Up", mesh),
    pMesh,
    model.Up()/beta
);

// Shear wave speed
pointScalarField Us
(
    IOobject("Us", mesh),
    pMesh,
    model.Us()*beta
);

// Stretch
pointScalarField stretch = mech.stretch();

// Wave speed for time increment
pointScalarField Up_time
(
    IOobject("Up_time", mesh),
    pMesh,
    model.Up()
);
Up_time = Up_time/stretch;


Info << "Creating fields for gradient ..." << endl;

// Gradient class
gradientSchemes grad(mesh, Sf, V);

// Gradient of nodal linear momentum
tensorField lmGrad = grad.gradient(lm);

// Reconstruction of nodal linear momentum
vectorField lmM(mesh.edges().size());
vectorField lmP(mesh.edges().size());
grad.reconstruct(lm, lmGrad, lmM, lmP);


Info << "Creating fields for fluxes ..." << endl;

// Edge linear momentum
vectorList lmC = inter.pointToEdge(lm);

// Edge traction
vectorField tC = (inter.pointToEdge(P) & N) + (0.5*mech.Smatrix(Up, Us) & (lmP-lmM));

Info << "Creating fields for angular momentum ..." << endl;

// Angular momentum class
angularMomentum am(mesh, mechanicalProperties);

// RHS of linear momentum equation
pointVectorField rhsLm
(
    IOobject ("rhsLm", mesh),
    pMesh,
    dimensionedVector("rhsLm", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
);

// RHS of linear momentum equation (stage 1)
pointVectorField rhsLm1 = rhsLm;

// RHS of angular momentum equation
pointVectorField rhsAm
(
    IOobject ("rhsAm", mesh),
    pMesh,
    dimensionedVector("rhsAm", dimensionSet(1,-1,-2,0,0,0,0), vector::zero)
);


Info << "Creating fields for post-processing ..." << endl;

// Nodal displacements
pointVectorField u
(
    IOobject ("u", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    pMesh,
    dimensionedVector("u", dimLength, vector::zero)
);

u.write();
X -= u;


Info << "Creating variables for time ..." << endl;

scalar tstep = 0;
dimensionedScalar deltaT("deltaT", dimTime, runTime.deltaTValue());
dimensionedScalar t("t", dimTime, 0.0);

if (timeStepping == "variable")
{
    deltaT = (cfl*h)/model.Up();
    runTime.setDeltaT(deltaT);
}

scalarList RKstage(2);
RKstage[0] = 0;
RKstage[1] = 1;


Info << "Printing data ..." << endl;

// Print material properties
model.printMaterialProperties();

// // Print global linear and angular momentum
am.printGlobalMomentum(lm,x);

// Print centroid of mesh
dual.printPrimalMeshCentroid();
dual.printDualMeshCentroid();