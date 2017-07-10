clear;clc
%% 
infile  = 'Data\precip.mon.ltm.v401.nc';
data    = nc_varget(infile, 'precip');
dataone = squeeze(data(1,:,:));
 

guess = 1;
gtype = 1;
nscan = 2000;
epsx  = 0.01;
relc  = 0.6;

tic
[datax,mscan,resmax] = fill_msg_grid(dataone, guess, gtype, nscan, epsx, relc);
toc