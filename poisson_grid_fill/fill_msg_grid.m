%########################################################################
% This is the entry subroutine
%########################################################################
function [out,mscan,resmax]= fill_msg_grid ( xio, guess, gtype, nscan, epsx, relc)     
   %integer             nscan, guess, gtype, mscan
   %double precision    xio(mx,ny), xmsg, epsx, relc

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % xout = fill_msg_grid(xio, gtype, nscan, eps, guess, opt)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   %if no missing values return the original array
   out    = xio;
   mscan  = 0;
   resmax = 0.0;
  
   miss = isnan(xio);
   nmsg = sum(miss(:));	  
   if (nmsg==0) 
       return 
   end

   [out,mscan,resmax] = poisxy2(xio,nscan,epsx,relc,guess,gtype);

   return 


%########################################################################
% This is real interpolation subroutine
%########################################################################
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

    %integer    maxscn, guess, gtype, mscan
    %double precision   a(il,jl), amsg, crit, resmax

    %%%%%%%%%% local %%%%%%%%%%%%%%%%
    %sor     = scratch area

    p25     = 0.25;
    [il,jl] = size(a);
    sor     = zeros(size(a));

    for j = 1 : jl
        n    = 0;
        aavg = 0.0d0;
        for i = 1 : il
            if (isnan(a(i,j)))
                sor(i,j) = relc;
            else
                n        = n + 1;
                aavg     = aavg + a(i,j);
                sor(i,j) = 0.0;
            end
        end

        if (n>0) 
            aavg = aavg/n;
        end
    
        if (guess==0) 
            for i = 1 : il
                if (isnan(a(i,j))) 
                    a(i,j) = 0.0;  
                end
            end
        elseif (guess==1) 
            for i = 1 : il
                if ( isnan(a(i,j))) 
                    a(i,j) = aavg;   
                end
            end 
        end 

	end % End of for J

  %-----------------------------------------------------------------------
  %     iterate until errors are acceptable.
  %-----------------------------------------------------------------------
	mscan = 0;  
    
	for it = 1 : maxscn
        resmax = 0.0;       
        mscan  = mscan + 1;
        for j = 1 : jl
            jp1 = j+1;
            jm1 = j-1;
        	
            if (j==1) 
                jm1 = 2;
            end
            
            if (j==jl)
            	jp1 = jl-1;
            end

            for i = 1 : il
                % only work on missing value locations
                if (sor(i,j)~=0.0) 
                    im1 = i-1;
                    ip1 = i+1;

                    % cyclic in x           
                    if (i==1  && gtype==1) 
                        im1 = il ;
                    end
                    if (i==il && gtype==1) 
                        ip1 = 1;
                    end

                    % not cyclic in x           
                    if (i==1  && gtype==0) 
                        im1 = 2;
                    end
                    if (i==il && gtype==0) 
                        ip1 = il-1; 
                    end

                    res    = p25*(a(im1,j)+a(ip1,j)+a(i,jm1)+a(i,jp1)) -a(i,j);
                    res    = res*sor(i,j);
                    a(i,j) = a(i,j) + res;
                    resmax  = max(abs(res),resmax);
                end
            end            

        end 
        
        out = a;
        
        % satisfy the condition, just stop
        if(resmax<= crit)           
            fprintf( 'Extrapolated  mscan= %4.4i scans and max residual= %14.7f\n',  mscan, resmax)
            break
        end
    end
    	
return
