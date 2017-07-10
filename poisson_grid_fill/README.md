# poisson_grid_fill

This method is forked from NCL. The original codes in formation F77 are translated into F90. 
In addition, add a MATLAB version. 

Description:

Replaces all missing values in a grid with values derived from solving Poisson's equation via relaxation.
The values at non-missing locations are used as boundary conditions and are returned unchanged.

Prototype - F90

    subroutine poisxy2(a,il,jl,amsg,maxscn,crit,relc,guess,gtype,resmax, mscan)
    !==================================================================================
    ! inputs:
    !     a       = array with missing areas to be filled with size of il*jl. 
    !     il      = number of points along 1st dimension to be filled
    !     jl      = number of points along 2nd dimension to be filled
    !     amsg    = missing value
    !     maxscn  = maximum number of passes allowed in relaxation
    !     crit    = criterion for ending relaxation before "maxscn" limit
    !     relc    = relaxation constant. Usually, 0.45 <= relc <= 0.6.
    !     gtype   = 0 : not cyclic in x
    !               1 : cyclic in x
    !     guess   = 0 : use 0.0 as an initial guess
    !             = 1 : at each "y" use the average values for that "y"
    !                   think zonal averages
    ! outputs:
    !
    !     a       = array with interpolated values 
    !               non missing areas remain unchanged.
    !     resmax  = max residual
    !     mscan   = real iteration times
    !==================================================================================


Prototype - MATLAB

		function [out,mscan,resmax]= poisxy2(a,maxscn,crit,relc,guess,gtype)
		    %%%%%%%%% inputs %%%%%%%%:
		    %a       = array with missing areas to be filled. 
		    %il      = number of points along 1st dimension to be filled
		    %jl      = number of points along 2nd dimension to be filled
		    %maxscn  = maximum number of passes allowed in relaxation
		    %crit    = criterion for ending relaxation before "maxscn" limit
		    %relc    = relaxation constant
		    %gtype   = 0 : not cyclic in x
		    %          1 : cyclic in x
		    %guess   = 0 : use 0.0 as an initial guess
		    %        = 1 : at each "y" use the average values for that "y"
		    %              think zonal averages
		
		    %%%%%%%%% outputs %%%%%%%%:
		    %out     = array with interpolated values 
		    %          non missing areas remain unchanged.
		    %resmax  = max residual


Examples - MATLAB

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
