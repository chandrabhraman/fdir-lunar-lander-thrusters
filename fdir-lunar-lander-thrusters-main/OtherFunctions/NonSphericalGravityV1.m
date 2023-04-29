function NSG=NonSphericalGravityV1(PosInertial,DCM,S,C,maxdeg,Re,GM)
% - Modified from gravitysphericalharmonic
% Input - PosInertial in Kms in Inertial Frame. 
         %DCM Converting ECI2ECEF
         %S-CoeffNormalized
         %C-CoeffNOmralized
         %maxdeg- Degree & Order
         %Re-Eq. Radius - Km.
         %GM - Gravitational Parameter - Km^3/Sec^2
      

Pos= DCM*PosInertial;   % Position in Fixed Frame. 

NormPos=norm(Pos);      % Geocentric Radius.

phic= asin(Pos(3)/NormPos); % Compute Geocentric Latitude.

lambda = atan2(Pos(2),Pos(1));

smlambda=zeros(maxdeg);
cmlambda=zeros(maxdeg);

slambda = sin(lambda);
clambda = cos(lambda);
smlambda(1) = 0;
cmlambda(1) = 1;
smlambda(2) = slambda;
cmlambda(2) = clambda;

for m=3:maxdeg+1
    smlambda(m) = 2.0.*clambda.*smlambda(m-1) - smlambda(m-2);
    cmlambda(m) = 2.0.*clambda.*cmlambda(m-1) - cmlambda(m-2);
end


% Compute normalized associated legendre polynomials
P = zeros(maxdeg+3, maxdeg+3);
scaleFactor = zeros(maxdeg+3, maxdeg+3);
cphi = cos(pi/2-phic);
sphi = sin(pi/2-phic);

% Seeds for recursion formula
P(1,1) = 1;            % n = 0, m = 0;
P(2,1) = sqrt(3)*cphi; % n = 1, m = 0;
scaleFactor(1,1) = 0;
scaleFactor(2,1) = 1;
P(2,2) = sqrt(3)*sphi; % n = 1, m = 1;
scaleFactor(2,2) = 0;

for n = 2:maxdeg+2
    k = n + 1;
    for m = 0:n
        p = m + 1;
        % Compute normalized associated legendre polynomials, P, via recursion relations 
        % Scale Factor needed for normalization of dUdphi partial derivative
                
        if (n == m)           
            P(k,k) = sqrt(2*n+1)/sqrt(2*n)*sphi*P(k-1,k-1);
            scaleFactor(k,k) = 0;
        elseif (m == 0)
            P(k,p) = (sqrt(2*n+1)/n)*(sqrt(2*n-1)*cphi*P(k-1,p) - (n-1)/sqrt(2*n-3)*P(k-2,p));
            scaleFactor(k,p) = sqrt( (n+1)*(n)/2);
        else
            P(k,p) = sqrt(2*n+1)/(sqrt(n+m)*sqrt(n-m))*(sqrt(2*n-1)*cphi*P(k-1,p) - sqrt(n+m-1)*sqrt(n-m-1)/sqrt(2*n-3)*P(k-2,p));
            scaleFactor(k,p) = sqrt( (n+m+1)*(n-m));
        end
    end
end

% Compute gravity in ECEF coordinates
[gx,gy,gz] = loc_gravityPCPF( Pos, maxdeg, P, C( 1:maxdeg+1, 1:maxdeg+1 ),S( 1:maxdeg+1, 1:maxdeg+1 ), smlambda,cmlambda, GM, Re, NormPos,scaleFactor );

NSG=DCM'*[gx;gy;gz];
end

function [gx,gy,gz] = loc_gravityPCPF(p,maxdeg,P,C,S,smlambda,cmlambda,GM,Re,r,scaleFactor)
rRatio   = Re/r;
rRatio_n = rRatio;

% initialize summation of gravity in radial coordinates
dUdrSumN      = 1;
dUdphiSumN    = 0;
dUdlambdaSumN = 0;

% summation of gravity in radial coordinates
for n = 2:maxdeg
    k = n+1;
    rRatio_n      = rRatio_n*rRatio;
    dUdrSumM      = 0;
    dUdphiSumM    = 0;
    dUdlambdaSumM = 0;
    for m = 0:n
        j = m+1;
        dUdrSumM      = dUdrSumM + P(k,j)*(C(k,j)*cmlambda(j) + S(k,j)*smlambda(j)); 
        dUdphiSumM    = dUdphiSumM + ( (P(k,j+1)*scaleFactor(k,j)) - p(3)/(sqrt(p(1)^2 + p(2)^2))*m*P(k,j))*(C(k,j)*cmlambda(j) + S(k,j)*smlambda(j)); 
        dUdlambdaSumM = dUdlambdaSumM + m*P(k,j)*(S(k,j)*cmlambda(j) - C(k,j)*smlambda(j));
    end
    dUdrSumN      = dUdrSumN      + dUdrSumM*rRatio_n*k;
    dUdphiSumN    = dUdphiSumN    + dUdphiSumM*rRatio_n;
    dUdlambdaSumN = dUdlambdaSumN + dUdlambdaSumM*rRatio_n;
end
% gravity in spherical coordinates
dUdr      = -GM/(r*r)*dUdrSumN;
dUdphi    =  GM/r*dUdphiSumN;
dUdlambda =  GM/r*dUdlambdaSumN;

% gravity in ECEF coordinates
gx = ((1/r)*dUdr - (p(3)/(r*r*sqrt(p(1)^2 + p(2)^2)))*dUdphi)*p(1)- (dUdlambda/(p(1)^2 + p(2)^2))*p(2); 
gy = ((1/r)*dUdr - (p(3)/(r*r*sqrt(p(1)^2 + p(2)^2)))*dUdphi)*p(2)+ (dUdlambda/(p(1)^2 + p(2)^2))*p(1); 
gz = (1/r)*dUdr*p(3) + ((sqrt(p(1)^2 + p(2)^2))/(r*r))*dUdphi;

% Special case for poles
atPole = abs(atan2(p(3),sqrt(p(1)^2 + p(2)^2)))==pi/2;
if any(atPole)
    gx(atPole) = 0;
    gy(atPole) = 0;
    gz(atPole) = (1/r(atPole))*dUdr(atPole)*p((atPole),3);
end

end
