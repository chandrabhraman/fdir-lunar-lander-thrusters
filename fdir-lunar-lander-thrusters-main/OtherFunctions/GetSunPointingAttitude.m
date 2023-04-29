% Get Sun Pointing Attitude 
function Q=GetSunPointingAttitude(Pos,Vel,JD,CentralBody)

% Input - Position & Velocity in Central Body J2000Frame - 3x1 Format
% JD - Julian Date
% CentralBody Name 'Eart' or 'Moon'

% Output = Quaternions - Aligned & Constraint - 
%          Aligned Solar Panel  '-XAxis' to Sun.
%          Constraint Vector - 'YAxis' to Orbit Normal

% Get Position of Sun from Central Body in J2000 Frame
et=cspice_str2et(sprintf('JD%0.9f',JD));
if(isrow(Pos))
    Pos=Pos';
end
if(isrow(Vel))
    Vel=Vel';
end
if(strcmp(CentralBody,'Eart'))
    PSun=cspice_spkpos('sun',et,'J2000','NONE','earth');
    if(isrow(PSun))
        PSun=PSun';
    end    
    PSun=PSun/norm(PSun);
elseif(strcmp(CentralBody,'Moon'))
    PSun=cspice_spkpos('sun',et,'J2000','NONE','moon');
    if(isrow(PSun))
        PSun=PSun';
    end    
    PSun=PSun/norm(PSun);
else
    disp('Error-CentralBodyNotFound')
end

P=Pos/norm(Pos)';   % Position Unit Vector
V=Vel/norm(Vel)';   % Velocity Unit Vector
N=cross(P,V)/sin(acos(dot(P,V)));       % Orbit Normal Unit Vector

SCap = PSun;
% MCap= cross(PSun,N)/norm(cross(PSun,N));
MCap= cross(PSun,N)/sin(acos(dot(PSun,N))); 
SSCap=[-0.965925826289068;0;0.258819045102521]; % Sun Vector in Body Frame 15 Deg from Negative X Axis
% SSCap=[-1;0;0];   % Sun Vector in Body Frame Negative X Axis 
R2=[0;1;0]; % Orbit Normal in Body Frame
MMCap = cross(SSCap,R2)/sin(acos(dot(SSCap,R2))); 

A= [SCap MCap cross(SCap,MCap)]*[SSCap MMCap cross(SSCap,MMCap)]';

Q=qGetQModified(A);
end

% Refer Triad Algorithm 
% https://en.wikipedia.org/wiki/Triad_method