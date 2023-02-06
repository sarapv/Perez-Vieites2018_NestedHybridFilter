function d = xyz2d ( m, x, y, z )

%*****************************************************************************80
%
%% XY2D converts a 2D Cartesian coordinate to a 1D Hilbert coordinate.
%
%  Discussion:
%
%    It is assumed that a square has been divided into an NxN array of cells,
%    where N=2^M is a power of 2.
%
%    Cell (0,0) is in the lower left corner, and (N-1,N-1) in the upper 
%    right corner.
%
%  Modified:
%
%    05 October 2017
%
%  Parameters:
%
%    Input, integer M, the index of the Hilbert curve.
%    The number of cells is N=2^M.
%    0 < M.
%
%    Input, integer X, Y, the Cartesian coordinates of a cell.
%    0 <= X, Y < N.
%
%    Output, integer D, the Hilbert coordinate of the cell.
%    0 <= D < N * N.
%
  if ( m <= 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'XY2D - Fatal error!\n' );
    fprintf ( 1, '  0 < M required.\n' );
    fprintf ( 1, '  M = %d\n', m );
    error ( 'XYZ2D: Fatal error!' );
  end

  n = m^3;

  if ( x < 0 || n <= x )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'XY2D - Fatal error!\n' );
    fprintf ( 1, '  0 <= X < N required.\n' );
    fprintf ( 1, '  N = %d\n', n );
    fprintf ( 1, '  X = %d\n', x );
    error ( 'XYZ2D: Fatal error!' );
  end

  if ( y < 0 || n <= y )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'XY2D - Fatal error!\n' );
    fprintf ( 1, '  0 <= Y < N required.\n' );
    fprintf ( 1, '  N = %d\n', n );
    fprintf ( 1, '  Y = %d\n', y );
    error ( 'XYZ2D: Fatal error!' );
  end

  if ( z < 0 || n <= z )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'XY2D - Fatal error!\n' );
    fprintf ( 1, '  0 <= Z < N required.\n' );
    fprintf ( 1, '  N = %d\n', n );
    fprintf ( 1, '  Z = %d\n', y );
    error ( 'XYZ2D: Fatal error!' );
  end
  
  xcopy = x;
  ycopy = y;
  zcopy = z;
  
  d = 0;

  s = floor ( m / 2 );

while ( 0 < s )

    if ( 0 < bitand ( uint32 ( abs ( xcopy ) ), uint32 ( s ) ) )
      rx = 1;
    else
      rx = 0;
    end

    if ( 0 < bitand ( uint32 ( abs ( ycopy ) ), uint32 ( s ) ) )
      ry = 1;
    else
      ry = 0;
    end

    if ( 0 < bitand ( uint32 ( abs ( zcopy ) ), uint32 ( s ) ) )
      rz = 1;
    else
      rz = 0;
    end
    
    d = d + s * s * s * ( bitxor( bitxor ( uint32 ( 7 * rx ), uint32 ( 3 * ry ) )  , uint32(rz)  )  );
    
    [ xcopy, ycopy, zcopy ] = rot3d ( s, xcopy, ycopy, zcopy, rx, ry, rz );

    s = floor ( s / 2 );

  end

  return
end