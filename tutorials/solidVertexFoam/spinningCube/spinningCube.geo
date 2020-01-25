//INPUTS

c = 20;

// Number of cells
cx = c;
cy = c;
cz = c;

// Lengths
lx = 1;
ly = 1;
lz = cz/cx;


Point(1) = {-lx/2, -ly/2, -lz/2};

Extrude {lx,0,0} {
     Point{1}; Layers{cx}; //Recombine;
}

Extrude {0,ly,0} {
     Line{1}; Layers{cy}; //Recombine;
}

Extrude {0,0,lz} {
     Surface{5}; Layers{cz}; //Recombine;
}


// Definition of surfaces for boundary conditions
Physical Surface("free") = {5,27,26,18,22,14};

//Physical Surface("free") = {26,18,22,14};
//Physical Surface("symmetricZ") = {5,27};

// Definition of a volume
Physical Volume("volume") = {1};



/*
// OUTPUTS

// Definition of points
Point(1) = {xmin, ymin, zmin};
Point(2) = {xmax, ymin, zmin};
Point(3) = {xmax, ymax, zmin};
Point(4) = {xmin, ymax, zmin};

// Definition of lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// Definition of line loop and surface
Line Loop(1) = {1,2,3,4} ;
Plane Surface(1) = {1};

// Definition of number of nodes on lines
Transfinite Line{1,3} = nx;
Transfinite Line{2,4} = ny;

// Definition of the corners of the mesh and creation of a structured hexahedral mesh
//Transfinite Surface{1} = {1,2,3,4};
//Recombine Surface{1};

// Extrusion in the third dimension
Extrude{0,0,zmax-zmin}{
Surface{1};
Layers{nz-1};//Recombine;
}

// Definition of surfaces for boundary conditions
Physical Surface("free") 		= {1,13,17,21,25,26};

// Definition of a volume
Physical Volume("volume") = {1};
*/
