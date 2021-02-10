
% compute a compund POD from several trajectories which already underwent a POD (because they had so many timesteps)

x_28 = load('28_single/VCS_fields_TT_pod_x.txt');
x_30 = load('30_single/VCS_fields_TT_pod_x.txt');
%x_32 = load('32_single/VCS_fields_TT_pod_x.txt');

y_28 = load('28_single/VCS_fields_TT_pod_y.txt');
y_30 = load('30_single/VCS_fields_TT_pod_y.txt');
%y_32 = load('32_single/VCS_fields_TT_pod_y.txt');

all_x = [x_30' x_28'];

%[U,S,V] = svd(x_30);

[U,S,V] = svd(all_x);
dim_x = size(all_x, 2);

Ut = U';

writematrix(Ut(1:dim_x,:), 'all_x.txt', 'Delimiter', ' ')
writematrix(dim_x);

diag(S(1:dim_x,1:dim_x))'

%writematrix(U(:,1:dim_x), 'all_x.txt', 'Delimiter', ' ')

all_y = [y_30' y_28'];

[U,S,V] = svd(all_y);
dim_y = size(all_y, 2);

diag(S(1:dim_y,1:dim_y))'

Ut = U';

writematrix(Ut(1:dim_y,:), 'all_y.txt', 'Delimiter', ' ')
writematrix(dim_y);


