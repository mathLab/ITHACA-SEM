x_30 = load('30_single/VCS_fields_TT_pod_x.txt');
x_32 = load('32_single/VCS_fields_TT_pod_x.txt');

y_30 = load('30_single/VCS_fields_TT_pod_y.txt');
y_32 = load('32_single/VCS_fields_TT_pod_y.txt');

all_x = [x_30' x_32'];

%[U,S,V] = svd(x_30);

[U,S,V] = svd(all_x);
dim_x = size(all_x, 2);
writematrix(U(:,1:dim_x), 'all_x.txt', 'Delimiter', ' ')

diag(S(1:45,1:45))'

all_y = [y_30' y_32'];

[U,S,V] = svd(all_y);
dim_y = size(all_y, 2);

diag(S(1:47,1:47))'

writematrix(U(:,1:dim_y), 'all_y.txt', 'Delimiter', ' ')



