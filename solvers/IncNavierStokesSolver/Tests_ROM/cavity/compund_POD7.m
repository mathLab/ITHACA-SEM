
% compute a compund POD from several trajectories which already underwent a POD (because they had so many time steps)

 x_100 = load('100_triple/VCS_fields_TT_pod_x.txt');
 x_110 = load('110_triple/VCS_fields_TT_pod_x.txt');
x_120 = load('120_triple/VCS_fields_TT_pod_x.txt');
% x_130 = load('130_triple/VCS_fields_TT_pod_x.txt');
x_140 = load('140_triple/VCS_fields_TT_pod_x.txt');
 x_150 = load('150_triple/VCS_fields_TT_pod_x.txt');

y_100 = load('100_triple/VCS_fields_TT_pod_y.txt');
y_110 = load('110_triple/VCS_fields_TT_pod_y.txt');
y_120 = load('120_triple/VCS_fields_TT_pod_y.txt');
% y_130 = load('130_triple/VCS_fields_TT_pod_y.txt');
y_140 = load('140_triple/VCS_fields_TT_pod_y.txt');
y_150 = load('150_triple/VCS_fields_TT_pod_y.txt');

all_x = [x_100' x_110' x_120' x_140' x_150'];
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


all_y = [y_100' y_110' y_120' y_140' y_150'];
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


