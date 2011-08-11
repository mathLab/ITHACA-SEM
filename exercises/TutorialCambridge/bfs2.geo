lcStep  = 0.1;
lcStep2 = 0.1;
lcInflow= 0.5;
lcFar   = 0.4;

Point(1) = {0,0,0,lcStep};
Point(2) = {0,-1,0,lcInflow};
Point(3) = {40,-1,0,lcFar};
Point(4) = {40,1,0,lcFar};
Point(5) = {-2,1,0,lcInflow};
Point(6) = {-2,0,0,lcInflow};

Point(10) = {-0.3,0,0,lcInflow};
Point(11) = {-0.3,1,0,lcInflow};
Point(12) = {5,-1,0,0.5};
Point(13) = {5,1,0,0.5};

Line(110) = {1,2};
Line(111) = {2,12};
Line(112) = {12,13};
Line(113) = {13,11};
Line(114) = {11,10};
Line(115) = {10,1};

Line(116) = {12,3};
Line(117) = {3,4};
Line(118) = {4,13};
Line(119) = {11,5};
Line(120) = {5,6};
Line(121) = {6,10};

Line Loop(220) = {110,111,112,113,114,115};
Line Loop(221) = {116,117,118,-112};
Line Loop(222) = {119,120,121,-114};

Transfinite Line{116,-118} = 45 Using Progression 1;
Transfinite Line{112,117,120,114} = 5 Using Bump 0.1;
Transfinite Line{119,-121} = 8 Using Progression 1.2;

Field[1] = Box;
Field[1].VIn = lcStep2;
Field[1].VOut = lcFar;
Field[1].XMax = 2;
Field[1].XMin = -0.2;
Field[1].YMin = -0.3;
Field[1].YMax = 0.3;
Field[1].ZMin = 0;
Field[1].ZMax = 0;



Background Field = 1;

Plane Surface(330) = {220};
Plane Surface(331) = {221};
Plane Surface(332) = {222};
Transfinite Surface{331} = {12,3,4,13};
Transfinite Surface{332} = {11,5,6,10};

//Recombine Surface{330};
Recombine Surface{331};
Recombine Surface{332};
Mesh.SubdivisionAlgorithm=2;
//Mesh.Triangles=1;

Transfinite Line {121, 119} = 8 Using Progression 1;
Transfinite Line {114, 120} = 5 Using Bump 0.1;

Physical Surface(333) = {332, 330, 331};
Physical Line(100) = {119, 113, 118, 121, 115, 110, 111, 116};
Physical Line(102) = {120};
Physical Line(103) = {117};