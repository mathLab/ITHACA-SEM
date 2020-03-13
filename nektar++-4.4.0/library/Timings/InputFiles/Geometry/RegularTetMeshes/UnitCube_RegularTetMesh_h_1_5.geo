Point(1) = {0,0,0,0.1};
l[] = Extrude {1,0,0} {
      Point{1}; Layers{5};
};
s[] = Extrude {0,1,0} {
      Line{l[1]}; Layers{5};
};
v[] = Extrude {0,0,1} {
      Surface{s[1]}; Layers{5};
};

MyVol=0;
Physical Volume (MyVol) = {v[1]};
Physical Surface(1) = {5};
Physical Surface(2) = {14};
Physical Surface(3) = {18};
Physical Surface(4) = {22};
Physical Surface(5) = {26};
Physical Surface(6) = {27};
