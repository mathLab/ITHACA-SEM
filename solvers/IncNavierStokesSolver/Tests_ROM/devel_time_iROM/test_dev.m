
% compute a compund POD from several trajectories which already underwent a POD (because they had so many time steps)


x_pod = load('test_param/VCS_fields_TT_pod_x.txt');
x_fom = load('test_param/VCS_fields_TT_x.txt');
x_fom_local = load('test_param/VCS_fields_TT_x_local.txt');

% measure the difference

diff12 = norm(x_fom(1,:) - x_fom(2,:)) / norm(x_fom(1,:))

diff45 = norm(x_fom(4,:) - x_fom(5,:)) / norm(x_fom(4,:))

y_pod = load('test_param/VCS_fields_TT_pod_y.txt');
y_fom = load('test_param/VCS_fields_TT_y.txt');
y_fom_local = load('test_param/VCS_fields_TT_y_local.txt');

all_x = [x_fom_local'];

[U,S,V] = svd(all_x);
dim_x = size(all_x, 2);

sing_vals = diag(S(1:dim_x,1:dim_x))';
modes_taken_x = sum(cumsum(sing_vals) ./ sum(sing_vals) < 1.1);

Ut = U';

writematrix(Ut(1:modes_taken_x,:), 'all_x.txt', 'Delimiter', ' ')
writematrix(modes_taken_x,  'dim_x.txt');

all_y = [y_fom_local'];

[U,S,V] = svd(all_y);
dim_y = size(all_y, 2);

sing_vals = diag(S(1:dim_y,1:dim_y))';
modes_taken_y = sum(cumsum(sing_vals) ./ sum(sing_vals) < 1.1);

Ut = U';

writematrix(Ut(1:modes_taken_y,:), 'all_y.txt', 'Delimiter', ' ')
writematrix(modes_taken_y,  'dim_y.txt');

modes_taken_x
modes_taken_y
