Point(1) = {0,0,0, 1.0};
Point(2) = {1,0,0, 1.0};
Point(3) = {1,1,0, 1.0};
Point(4) = {0,1,0, 1.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(5) = {1,2,3,4};

Plane Surface(10) = {5};
