clear all; close all; clc; %#ok<CLALL>

ntime = 21;
nsd   = 3;
nNo   = 125;
tmin  = 0.0;
tmax  = 0.001;

u = zeros(nNo,nsd,ntime);
for i=1:ntime
    for j=1:nsd
        fname = sprintf('csv/fb%d_%02d.csv',j,i);
        V = csvread(fname)';
        u(:,j,i) = reshape(V, nNo, 1);
    end
end

t = linspace(tmin, tmax, ntime);
fname = 'bforce.dat';
fid = fopen(fname,'w');
fprintf(fid,'%d   %d   %d\n', nsd, ntime, nNo);
for i=1:ntime
    fprintf(fid,'%.6f\n', t(i));
end
for a=1:nNo
    fprintf(fid,'%d\n',a);
    for i=1:ntime
        for j=1:nsd
            fprintf(fid,'%.9f ', u(a,j,i));
        end
        fprintf(fid,'\n');
    end
end
fclose(fid);


