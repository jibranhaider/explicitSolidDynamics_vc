// Number of nodes
c = 4;      // 625 nodes
//c = 8;      // 3969 nodes
//c = 16;     // 28033 cells

// Lengths
lx = 1;
ly = 6;
lz = 1;

// Number of cells
cx = c;
cy = 6*c;
cz = c;


Point(1) = {-lx/2, 0, -lz/2};

Extrude {lx,0,0} {
     Point{1}; Layers{cx};
}

Extrude {0,ly,0} {
     Line{1}; Layers{cy};
}

Extrude {0,0,lz} {
     Surface{5}; Layers{cz};
}

// Definition of surfaces for boundary conditions
Physical Surface("free") = {5,27,26,18,22};
Physical Surface("bottom") = {14};

// Definition of a volume
Physical Volume("volume") = {1};
