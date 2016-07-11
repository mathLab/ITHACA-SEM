// Quad mesh of a square
//
// Run: gmsh -2 -order n SquareQuad.geo

Point(1) = {0,0,0,0.1};
l[] = Extrude {1,0,0} {
      Point{1}; Layers{4};
};
s[] = Extrude {0,1,0} {
      Line{l[1]}; Layers{4};
};
Physical Surface(0) = {5};
Physical Line(1) = {1};
Physical Line(2) = {4};
Physical Line(3) = {2};
Physical Line(4) = {3};
Recombine Surface {5};
Recombine Surface {5};
Recombine Surface {5};
Recombine Surface {5};
