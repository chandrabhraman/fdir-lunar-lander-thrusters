% Source GMAT Math Spec
function [a,e,i,AOP,Omega,TA]=Cartesian2Keplerian(Position,Velocity,CentralBody)
% Input - Position - Km, Velocity - Km
% Input - CentralBody , 'Eart', 'Moon'

% Output 
% a - SemiMajorAxis -Km
% e - Eccentricity , i - Inclination - deg, AOP - deg, Omega - RAAN - deg,
% TA - deg.
if(isrow(Position))
    Position=Position';
end

if(isrow(Velocity))
    Velocity=Velocity';
end

if(strcmp(CentralBody,'Eart'))
    mu=398600.4415;
elseif(strcmp(CentralBody,'Moon'))
    mu=4902.8005821478;
    
end

% Specific Angular Momentum 
h=cross(Position,Velocity);
    hMag = norm(h);
    
% Vector in the Direction of Line of Nodes

n = cross([0;0;1],h);
    nMag = norm(n);
    
RMag=norm(Position);
VMag=norm(Velocity);
    EVector = ((VMag^2 - mu/RMag)*Position - dot(Position,Velocity)*Velocity)/mu;
    e = norm(EVector);  % Eccentricity
    
    Energy = VMag^2/2 - mu/RMag;
    
    if(abs(1-e)<10^-7)
        disp('Error - Orbit Near Parabolic')
    else
        a = -mu/(2*Energy);  % SemiMajor Axis
        if(abs(a*(1-e))<0.001)
            disp('Error - Conic Section is Nearly Singular')
        else
             i = acos(h(3)/hMag);    % Inclination
        end
    end
    % Non-circular, Inclined Orbit
    if(e>= 10^-11 && i>=10^-11 && i<=(pi-10^-11))
        Omega = acos(n(1)/nMag);
        if(n(3)<0)
            Omega = 2*pi - Omega;
        end
        AOP = acos(dot(n,EVector)/(nMag*e));
        if(EVector(3)<0)
            AOP = 2*pi - AOP;
        end
        TA = acos(dot(EVector,Position)/(e*RMag));
        if(dot(Position,Velocity)<0)
            TA = 2*pi - TA;
        end
    end
    % Non-circular, Equatorial Orbit
    if(e>= 10^-11 && (i<10^-11 || i>(pi-10^-11)))
        Omega = 0;
        AOP = acos(EVector(1)/e);
        if(EVector(2)<0)
            AOP = 2*pi - AOP;
        end
        TA = acos(dot(EVector,Position)/(e*RMag));
        if(dot(Position,Velocity)<0)
            TA = 2*pi -TA;
        end
    end
    % Circular, Inclined Orbit
     if(e< 10^-11 && i>=10^-11 && i<=(pi-10^-11))
         Omega = acos(n(1)/nMag);
         if(n(2)<0)
             Omega = 2*pi - Omega;
         end
         AOP =0;
         TA = acos(dot(n,Position)/(nMag*RMag));
         if(Position(3)<0)
             TA = 2*pi - TA;
         end
     end
   % Circular, Equatorial Orbit
    if(e< 10^-11 && (i<10^-11 || i>(pi-10^-11)))
        Omega =0;
        AOP = 0;
        TA = acos(Position(1)/RMag);
        if(Position(2)<0)
            TA = 2*pi - TA;
        end
    end
    
    TA = real(TA * 180/pi);
    AOP = real(AOP * 180/pi);
    Omega = real(Omega * 180/pi);
    i = real(i*180/pi);
    
end


    
   
        
    
   
  
   
   
   
   
   