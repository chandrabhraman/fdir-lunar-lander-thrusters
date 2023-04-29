function Q = qGetQModified( R )
% qGetQ: converts 3x3 rotation matrix into equivalent quaternion
% Q = qGetQ( R );

% NOTE:[27-11-14] Properties noted for this algorithm: One component of output q will
% definitely be positive. The positive sign is accorded to the q component
% having the max modulo. The scheme identifies individual component modulos
% by working with the diagonal entries and using the quaternion 2 norm =1
% property. When the sign for the max modulo component is assigned, the
% rest of the terms are assigned signs by manipulating off diagonal terms
% of the rotation matrix. Hence, the problem with using this scheme is:
% The sign of the quaternions might switch if there is a new element having
% a max modulo. 

% Also important is the fact that the rotation matrix corresponds not to
% the rotation matrix which converts a vector in the base frame to the
% current frame but the transpose of the same. 


%%%% Input - Roation Matrix from Body To Base Frame

%%%% Ouput - Quaternions with Q4-Scalar

[r,c] = size( R );
if( r ~= 3 | c ~= 3 )
    fprintf( 'R must be a 3x3 matrix\n\r' );
    return;
end

Rxx = R(1,1); Rxy = R(1,2); Rxz = R(1,3);
Ryx = R(2,1); Ryy = R(2,2); Ryz = R(2,3);
Rzx = R(3,1); Rzy = R(3,2); Rzz = R(3,3);

w = sqrt( trace( R ) + 1 ) / 2;

% check if w is real. Otherwise, zero it.
if( imag( w ) > 0 )
     w = 0;
end

x = sqrt( 1 + Rxx - Ryy - Rzz ) / 2;
y = sqrt( 1 + Ryy - Rxx - Rzz ) / 2;
z = sqrt( 1 + Rzz - Ryy - Rxx ) / 2;

[element, i ] = max( [w,x,y,z] );

global icount
icount=[icount i];

if( i == 1 )
    x = ( Rzy - Ryz ) / (4*w);
    y = ( Rxz - Rzx ) / (4*w);
    z = ( Ryx - Rxy ) / (4*w);
end

if( i == 2 )
    w = ( Rzy - Ryz ) / (4*x);
    y = ( Rxy + Ryx ) / (4*x);
    z = ( Rzx + Rxz ) / (4*x);
end

if( i == 3 )
    w = ( Rxz - Rzx ) / (4*y);
    x = ( Rxy + Ryx ) / (4*y);
    z = ( Ryz + Rzy ) / (4*y);
end

if( i == 4 )
    w = ( Ryx - Rxy ) / (4*z);
    x = ( Rzx + Rxz ) / (4*z);
    y = ( Ryz + Rzy ) / (4*z);
end



Q = [ x; y; z; w];