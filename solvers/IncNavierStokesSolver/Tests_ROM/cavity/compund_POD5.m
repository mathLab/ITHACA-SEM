
% compute a compund POD from several trajectories which already underwent a POD (because they had so many time steps)

%x_28 = load('28_single/VCS_fields_TT_pod_x.txt');
%x_30 = load('30_single/VCS_fields_TT_pod_x.txt');
%x_32 = load('32_single/VCS_fields_TT_pod_x.txt');
%x_130 = load('130_triple/VCS_fields_TT_pod_x.txt');
x_170 = load('170_triple/VCS_fields_TT_pod_x.txt');

%y_28 = load('28_single/VCS_fields_TT_pod_y.txt');
%y_30 = load('30_single/VCS_fields_TT_pod_y.txt');
%y_32 = load('32_single/VCS_fields_TT_pod_y.txt');
%y_130 = load('130_triple/VCS_fields_TT_pod_y.txt');
y_170 = load('170_triple/VCS_fields_TT_pod_y.txt');

%all_x = [x_30' x_28'];
all_x = [x_170'];

%[U,S,V] = svd(x_30);

[U,S,V] = svd(all_x);
dim_x = size(all_x, 2);

sing_vals = diag(S(1:dim_x,1:dim_x))';
modes_taken_x = sum(cumsum(sing_vals) ./ sum(sing_vals) < 0.99);

Ut = U';

writematrix(Ut(1:modes_taken_x,:), 'all_x.txt', 'Delimiter', ' ')
writematrix(modes_taken_x,  'dim_x.txt');

%diag(S(1:dim_x,1:dim_x))'

%writematrix(U(:,1:dim_x), 'all_x.txt', 'Delimiter', ' ')

%all_y = [y_30' y_28'];
all_y = [y_170'];

[U,S,V] = svd(all_y);
dim_y = size(all_y, 2);

%diag(S(1:dim_y,1:dim_y))'

sing_vals = diag(S(1:dim_y,1:dim_y))';
modes_taken_y = sum(cumsum(sing_vals) ./ sum(sing_vals) < 0.99);

Ut = U';

writematrix(Ut(1:modes_taken_y,:), 'all_y.txt', 'Delimiter', ' ')
writematrix(modes_taken_y,  'dim_y.txt');


