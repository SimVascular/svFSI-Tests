clear all; close all; clc; %#ok<CLALL>

ntime = 21;
nsd   = 3;
tmin  = 0.0;
tmax  = 0.001;

nFa   = 6;
flist = ['X0'; 'X1'; 'Y0'; 'Y1'; 'Z0'; 'Z1'];

for fa=1:nFa
    fhdr   = flist(fa,:);
    fname  = sprintf('csv/bc_%s_nodeid.csv', fhdr);
    nodeId = csvread(fname);
    nNo = size(nodeId,1);

    h = zeros(nNo,nsd,ntime);
    for i=1:ntime
        for j=1:nsd
            fname = sprintf('csv/bc_%s_h%d_%02d.csv',...
            fhdr,j,i);
            if ~exist(fname,'file')
                continue;
            end
            V = csvread(fname)';
            h(:,j,i) = reshape(V, nNo, 1);
        end
    end

    t = linspace(tmin, tmax, ntime);
    fname = sprintf('%s_hbc.dat', fhdr);
    fid = fopen(fname,'w');
    fprintf(fid,'%d   %d   %d\n', nsd, ntime, nNo);
    for i=1:ntime
        fprintf(fid,'%.6f\n', t(i));
    end
    for a=1:nNo
        fprintf(fid,'%d\n', nodeId(a));
        for i=1:ntime
            for j=1:nsd
                fprintf(fid,'%.9f ', h(a,j,i));
            end
            fprintf(fid,'\n');
        end
    end
    fclose(fid);
end

