lc = 0.5;
r_i = 1.0;
r_o = 1.2;

// Annular points, z = -1
Point(1)  = {0  , 0  , -1, lc}; /* Circle origins */
Point(2)  = {r_i, 0  , -1, lc};
Point(3)  = {0  , r_i, -1, lc};
Point(4)  = {r_o, 0  , -1, lc};
Point(5)  = {0  , r_o, -1, lc};

// Annular points, z = 1
Point(6)  = {0  , 0  ,  1, lc}; /* Circle origins */
Point(7)  = {r_i, 0  ,  1, lc};
Point(8)  = {0  , r_i,  1, lc};
Point(9)  = {r_o, 0  ,  1, lc};
Point(10) = {0  , r_o,  1, lc};

// Annulus, z = -1
Circle(1) = {2, 1, 3};
Line  (2) = {3, 5};
Circle(3) = {5, 1, 4};
Line  (4) = {4, 2};
Line Loop(1) = {1, 2, 3, 4};

// Annulus, z = 1
Circle(5) = {7, 6, 8};
Line  (6) = {8, 10};
Circle(7) = {10, 6, 9};
Line  (8) = {9, 7};
Line Loop(2) = {5, 6, 7, 8};

// Connecting lines
Line (10) =  {2, 7};
Line (11) =  {3, 8};
Line (12) =  {4, 9};
Line (13) =  {5, 10};

// End planar surfaces
Line Loop (5) = {11, 6, -13, -2};
Line Loop (6) = {10, -8, -12, 4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

// Ruled surfaces
Line Loop(3) = {1, 11, -5, -10};
Line Loop(4) = {3, 12, -7, -13};

Surface(1) = {1};
Surface(2) = {2};
Surface(3) = {3};
Surface(4) = {4};

// Surface loop for annulus extrusion.
Surface Loop(1) = {1,2,3,4,5,6};

// Define a bounding box to fit everything into.
xmin = -1.0;
xmax = r_o + 1.0;
ymin = -1.0;
ymax = r_o + 1.0;
zmin = -2.0;
zmax = 2.0;
Point(100) = {xmin, ymin, zmin, lc};
Point(101) = {xmax, ymin, zmin, lc};
Point(102) = {xmax, ymax, zmin, lc};
Point(103) = {xmin, ymax, zmin, lc};
Point(104) = {xmin, ymin, zmax, lc};
Point(105) = {xmax, ymin, zmax, lc};
Point(106) = {xmax, ymax, zmax, lc};
Point(107) = {xmin, ymax, zmax, lc};

Line (100) = {100, 101};
Line (101) = {101, 102};
Line (102) = {102, 103};
Line (103) = {103, 100};
Line (104) = {104, 105};
Line (105) = {105, 106};
Line (106) = {106, 107};
Line (107) = {107, 104};
Line (108) = {100, 104};
Line (109) = {101, 105};
Line (110) = {102, 106};
Line (111) = {103, 107};

Line Loop(100) = {100, 101,  102,  103}; // bottom
Line Loop(101) = {100, 109, -104, -108};
Line Loop(102) = {101, 110, -105, -109};
Line Loop(103) = {102, 111, -106, -110};
Line Loop(104) = {103, 108, -107, -111};
Line Loop(105) = {104, 105,  106,  107}; // top

Plane Surface(100) = {100};
Plane Surface(101) = {101};
Plane Surface(102) = {102};
Plane Surface(103) = {103};
Plane Surface(104) = {104};
Plane Surface(105) = {105};

Surface Loop(100) = { 100, 101, 102, 103, 104, 105 };

Volume(1) = {100, -1};
