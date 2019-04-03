clear all
close all
fid = fopen('./Diagr2/bif_diagr.txt');

fromfile=fscanf(fid,'%f %f');
fclose(fid);

dati=zeros(length(fromfile)/2,2);
j=1;
for i=1:2:length(fromfile)
    dati(j,[1,2])=fromfile([i,i+1]);
    j=j+1;
end


scatter(dati(:,1),dati(:,2),'filled');
xlabel("\nu");
ylabel("v(x_0,y_0)");