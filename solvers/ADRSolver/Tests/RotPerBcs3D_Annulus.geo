Point(1) = {1, 0, 0, 1.0};
Point(2) = {0.3, 0, 0, 1.0};
Line(3) = {1, 2};
Transfinite Line(3) = 3;
s[] = Extrude {0, 0, 1} {
  Line{3}; Layers{3}; Recombine;
};
v[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{s[1]}; Layers{10}; Recombine;
};

Physical Surface(1) = {24}; // inflow
Physical Surface(2) = {16}; // outflow
Physical Surface(3) = {28}; // outer wall
Physical Surface(4) = {20}; // inner wall
Physical Surface(5) = {7};  // face at 0 deg
Physical Surface(6) = {29}; // face at 90 deg
Physical Volume(0) = {1}; // volume
