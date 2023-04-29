% Position & Velocity Integration 
% http://www.nist.gov/pml/div688/grp50/leapsecond.cfm - Diff. b/w UT1-UTC
% http://www.leapsecond.com/java/gpsclock.htm - Diff. b/w TAI & UTC
% http://maia.usno.navy.mil/search/search.html - Polar Motion Values.
% https://hpiers.obspm.fr/eop-pc/   - Lastest Polar Motion, UTC1-UTC,
% TDB-UTC, 
function [PosVelDot]=FunPosVel(t,PosVelo,NSG,Other,ENV,FunPosVelInput)

JulianArray=Other.JulianArray; CentralBody=Other.CentralBody; mass_sc=Other.mass_sc; BodyAcc = Other.BodyAcc;
SineCoef =NSG.SineCoef; CosCoef = NSG.CosCoef; Re=NSG.Re; GM=NSG.GM; DeltaAT=NSG.DeltaAT; DeltaUTC=NSG.DeltaUTC; PolarX=NSG.PolarX; PolarY=NSG.PolarY; NutationData=NSG.NutationData; Degree=NSG.Degree; DeltaPsi=NSG.DeltaPsi; DeltaEpsilon=NSG.DeltaEpsilon;

JDT=JulianArray+t/86400;
Rx=PosVelo(1);
Ry=PosVelo(2);
Rz=PosVelo(3);
Vx=PosVelo(4);
Vy=PosVelo(5);
Vz=PosVelo(6);


PosInertial=[Rx;Ry;Rz];

et=cspice_str2et(sprintf('JD%0.9f',JDT)); %Get ephemeris Time from SPICE function

if(strcmp(CentralBody,'Eart'))
    DCMJ2000toFixed=FK5(JDT,DeltaAT,DeltaUTC,[DeltaPsi DeltaEpsilon],[PolarX PolarY],NutationData); % ECI2ECEF
    [psun]=cspice_spkpos('sun',et,'J2000','NONE','earth');
    [pmoon]=cspice_spkpos('moon',et,'J2000','NONE','earth');
    PointMassAcc=PointMassGravityEarth(PosInertial,psun,pmoon);     % Gravity Due to Moon and Sun
    [SRP,~] = shadowfunction(PosInertial,psun,6378.1363,mass_sc);
elseif(strcmp(CentralBody,'Moon'))
    JDTDB = utc2tdb(JDT,ENV.jdateleap,ENV.leapsec);
    [Libration,~] = jplephemSimulink(JDTDB, 15, 0,ENV);
    phi = Libration(1);
    theta = Libration(2);
    psi = mod(Libration(3), 2.0 * pi);
    DCMJ2000toFixed=angle2dcm(phi,theta,psi,'ZXZ');    % MJ2000 to MoonPA
    [psun]=cspice_spkpos('sun',et,'J2000','NONE','moon');
    [pearth]=cspice_spkpos('earth',et,'J2000','NONE','moon');
    PointMassAcc=PointMassGravityMoon(PosInertial,psun,pearth);     % Gravity Due to Moon and Sun
    [SRP,~] = shadowfunction(PosInertial,psun,1738.2,mass_sc);
end

     NSG=NonSphericalGravityV1(PosInertial,DCMJ2000toFixed,SineCoef,CosCoef,Degree,Re,GM);   % NonSpherical Gravity 
     
    
     if(FunPosVelInput.TargeterEnable==1)
          RECI  = [Rx Ry Rz];
          VECI = [Vx Vy Vz];
         if(FunPosVelInput.TCMEnable==0)
             if(strcmp(CentralBody,'Eart')==1)
                % Earth Burns -- ZCap Along Velocity Vector, YCap is Opp. to Orbit Normal, XCap Completes Right Hand Rule (VNC-Earth) 
                ZCap=VECI/norm(VECI);   % Unit Vector along Velocity Vector - For Earth Orbits - Body Z is along Velcoity Direction. 
                YCap=cross(ZCap,RECI/norm(RECI))/sin(acos(dot(ZCap,RECI/norm(RECI))));   % For Earth Orbit -Y is Along Orbit Normal
                XCap=cross(YCap,ZCap);  % Orbital Normal Cross Velocity Vector.
                R=[XCap' YCap' ZCap'];  % Making the Rotation Matrix
                qb=qGetQModified(R);   % Getting Attitude(Quaternions) from Rotation Matrix
                qb=qnorm(qb); % Normalizing Quaternions
            elseif(strcmp(CentralBody,'Moon')==1)
                %Moon Burns -- ZCap Along Anti Velocity, XCap is Opp. to Orbit Normal, YCap Completes Right Hand Rule (VNC-Moon) 
                ZCap=-VECI/norm(VECI);   % Unit Vector along Velocity Vector - For Earth Orbits - Body Z is along Velcoity Direction. 
                XCap=cross(RECI/norm(RECI),ZCap)/sin(acos(dot(RECI/norm(RECI),ZCap)));   % For Earth Orbit -Y is Along Orbit Normal
                YCap=cross(ZCap,XCap);  % Orbital Normal Cross Velocity Vector.
                R=[XCap' YCap' ZCap'];  % Making the Rotation Matrix
                qb=qGetQModified(R);   % Getting Attitude(Quaternions) from Rotation Matrix
                qb=qnorm(qb); % Normalizing Quaternions
             end
           QMatrix=[qb(1)^2-qb(2)^2-qb(3)^2+qb(4)^2 2*(qb(1)*qb(2)+qb(3)*qb(4)) 2*(qb(1)*qb(3)-qb(2)*qb(4));2*(qb(1)*qb(2)-qb(3)*qb(4)) -qb(1)^2+qb(2)^2-qb(3)^2+qb(4)^2 2*(qb(2)*qb(3)+qb(1)*qb(4));2*(qb(1)*qb(3)+qb(2)*qb(4)) 2*(qb(2)*qb(3)-qb(1)*qb(4)) -qb(1)^2-qb(2)^2+qb(3)^2+qb(4)^2]; % Inertial to Body
           ThrustInertial=QMatrix'*[0;0;FunPosVelInput.ConstantThrust];
           ElapsedTime = t - FunPosVelInput.RegisterBurnStartTime; % Seconds
           CurrentMass = mass_sc - FunPosVelInput.MassFlowRate*ElapsedTime;
           BodyAcc=ThrustInertial/(CurrentMass*1000);
         elseif(FunPosVelInput.TCMEnable==1)
             XCap=VECI/norm(VECI);
             RCap=RECI/norm(RECI);
             NCap=cross(RCap,XCap)/sin(acos(dot(RCap,XCap)));
             BCap=cross(XCap,NCap);
             R=[XCap' NCap' BCap'];  % Body to Inerital 
%              UnitVectorVNB=[XDire;YDire;ZDire]/norm([XDire;YDire;ZDire]);
             XDireNew = cosd(FunPosVelInput.InPlane)*cosd(FunPosVelInput.OutPlane);   YDireNew = cosd(FunPosVelInput.InPlane)*sind(FunPosVelInput.OutPlane);    ZDireNew = sind(FunPosVelInput.InPlane);   % Added Newley
             UnitVectorVNB = [XDireNew;YDireNew;ZDireNew];  % Added Newely 
             % Added Newely - InPlane, OutPlane to Global Variables
             ThrustInertial=R*UnitVectorVNB*FunPosVelInput.ConstantThrust; 
             ElapsedTime = t - FunPosVelInput.RegisterBurnStartTime; % Seconds
             CurrentMass = mass_sc - FunPosVelInput.MassFlowRate*ElapsedTime;
             BodyAcc=ThrustInertial/(CurrentMass*1000);
         end
     end
     
     TotalAcc = NSG+PointMassAcc+BodyAcc;
     

GGx=TotalAcc(1);
GGy=TotalAcc(2);
GGz=TotalAcc(3);
RxDot=Vx;
RyDot=Vy;
RzDot=Vz;
VxDot=GGx;
VyDot=GGy;
VzDot=GGz;


PosVelDot=[RxDot;RyDot;RzDot;VxDot;VyDot;VzDot];
end