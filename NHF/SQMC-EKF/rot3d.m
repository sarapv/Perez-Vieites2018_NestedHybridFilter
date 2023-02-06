function [ x, y , z] = rot3d ( s, x, y, z, rx, ry, rz ) 

%*****************************************************************************80
%
%% ROT rotates and flips a quadrant appropriately.
%
%  Modified:
%
%    05 December 2015
%
%  Parameters:
%
%    Input, integer N, the length of a side of the square.  
%    N must be a power of 2.
%
%    Input/output, integer X, Y, the coordinates of a point.
%
%    Input, integer RX, RY, ???
%
    if ( ry == 0 )
%
%  Reflect.
%
        if ( rx == 1)
 
            if(rz == 1) % rx = 1, ry = 0, rz = 1
                  xaux = x-s;
                  
                  zaux = z-s;
                  xaux = bitxor(uint32(xaux),s-1);
                  yaux= bitxor(uint32(y),s-1);
                  zaux = bitxor(uint32(zaux),0);
                  t = xaux;
                  x = yaux;
                  y = zaux;
                  z = t;
                  
            else % rx = 1, ry = 0, rz = 0
                xaux = x-s;
                
                xaux = bitxor(uint32(xaux),s-1);
                yaux= bitxor(uint32(y),0);
                zaux = bitxor(uint32(z),s-1);
                t = zaux;
                y = xaux;
                z = yaux;
                x = t;
                
            end
        
        else 
        
%  Flip.
    

            if (rz == 1) % rx = ry = 0; rz = 1
                zaux = z - s;
                t =x;
                x = y;
                y = zaux;
                z = t;
            else    % rx = ry = rz = 0
                t = x;
                x = z;
                z = y;
                y = t;
            end
    

        end
    
    
    else
          
        if (rz == 0) % rx = 0 & 1 , ry = 1, rz = 0
            yaux = y - s;
            z = bitxor(z,s-1);
            y = bitxor(yaux,s-1);
        else
            
            if(rx == 0) %rx =0, ry = 1, rz = 1
                yaux = y - s;
                zaux = z - s;
                t =x;
                x = yaux;
                y = zaux;
                z = t;
            else  %rx =1, ry = 1, rz = 1
                xaux = x - s;
                yaux = y-s;
                zaux = z-s;
                xaux = bitxor(uint32(xaux),s-1);
                yaux= bitxor(uint32(yaux),s-1);
                zaux = bitxor(uint32(zaux),0);
                t = xaux;
                x = yaux;
                y = zaux;
                z = t;
            end
            
        end
        
        
          
    end

  return
end