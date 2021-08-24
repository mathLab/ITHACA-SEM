

max_reduc_int = 55;

all_data = [];

for i = 1:max_reduc_int
    filename = strcat('ROM_cluster_reduc', int2str(i) ,'.txt');
    dataset = load(filename);
    all_data = [all_data dataset];
end

plot(log10(mean(all_data)))

