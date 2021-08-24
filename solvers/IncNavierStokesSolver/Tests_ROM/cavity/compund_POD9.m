
% compute a compund POD from several trajectories which already underwent a POD (because they had so many time steps)

x_500 = load('500_triple/VCS_fields_TT_pod_x.txt');
x_700 = load('700_triple/VCS_fields_TT_pod_x.txt');

y_500 = load('500_triple/VCS_fields_TT_pod_y.txt');
y_700 = load('700_triple/VCS_fields_TT_pod_y.txt');

all_x = [ x_500' x_700'];
%all_x = [ x_120'  x_140' ];


[U,S,V] = svd(all_x);
dim_x = size(all_x, 2);

Ut = U';

sing_vals = diag(S(1:dim_x,1:dim_x))';
modes_taken_x = sum(cumsum(sing_vals) ./ sum(sing_vals) < 0.99);

modes_taken_x
dim_x

writematrix(Ut(1:modes_taken_x,:), 'all_x.txt', 'Delimiter', ' ')
writematrix(modes_taken_x, 'dim_x.txt');


all_y = [ y_500' y_700'];
%all_y = [ y_120'  y_140' ];

[U,S,V] = svd(all_y);
dim_y = size(all_y, 2);
sing_vals = diag(S(1:dim_y,1:dim_y))';
modes_taken_y = sum(cumsum(sing_vals) ./ sum(sing_vals) < 0.99);

modes_taken_y
dim_y
Ut = U';

writematrix(Ut(1:modes_taken_y,:), 'all_y.txt', 'Delimiter', ' ')
writematrix(modes_taken_y, 'dim_y.txt');


