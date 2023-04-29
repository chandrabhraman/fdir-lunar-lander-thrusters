function [Text,SolarTorque,MagTorque,GravTorque]=dist_torque(JD,t,qb,Isc,Planet,OrbitElements,rcg,SunVectorJ2000,APanel,BPanel,CPanel)

% Note: - Solar Radiation Torque Should be Changed. Center of Pressure
% Should be Determined from Sun Angle with Solar Panels.

%%%%%%%%

% Input:- Current Julian Date, Time Past Perigee, Attitude of Spacecraft in
% J2000 Frame, Current Inertia of Spacecraft, Central Body Name, Orbital
% Elements: Ecc, SMA, RAAN, AOP, Inclination,CurrentCGPositionInModelFrame.

% Output:- Disturbance Torque. 
% For Central Body Earth. Disturbance Torque- Gravitational Torque,
% Magnetic Torque, Atmospheric Torque, Solar Radiation Torque. 

% For Central Body Moon. Disturbace Torque- Gravitational Torque, Solar
% Radiation Torque.

%%%%%%%%

%%%%%%%%%%%
% Ver1 - 18-7-15 , Changed Magnetic Torque Model. Added Rotation's to Magnetic Frame 
% Ver2 - 23-7-15, Included Solar Panel Torques Reliastic Model
%%%%%%%%%%%%%%

if(strcmp(Planet,'Eart'))
    GConstant=398600;                % Earth - Km^3/sec^2
    PlanetRadius=6371;               % EarthRadius-Km
elseif(strcmp(Planet,'Moon'))
    GConstant=4902.8;                % Moon Gravitational Parameter-Km^3/sec^2
    PlanetRadius=1738;               % MoonRadius-Km
end
    ecc=OrbitElements(1);                         % Eccentricity
    SMA=OrbitElements(2);                         % Semi Major Axis.
if(SMA >0)    
    Timeperiod=2*pi*sqrt(SMA^3/GConstant);        % TimePeriod in Seconds.
    TAdot=degtorad((360/Timeperiod)*0.016);                                 % Rate of Change of TA for every 0.016 Seconds.
else
    TAdot = 0;  % For Hyperbolic Orbits TAdot cannot be calculated from above formula. So Temporary
end
    RAAN=degtorad(- OrbitElements(4));            % RANN in Radians.
    AOP=degtorad(-OrbitElements(5));              % AOP.
    i=degtorad(-OrbitElements(3));                % Inclination.

% Transforming RTN to ECI 
Alpha3=RAAN;
RotRANN=[cos(Alpha3) sin(Alpha3) 0;-sin(Alpha3) cos(Alpha3) 0;0 0 1];   % Rotation Around Z Axis-- RAAN.
Alpha1=i;
RotINC=[1 0 0;0 cos(Alpha1) sin(Alpha1);0 -sin(Alpha1) cos(Alpha1)];    % Rotation Around X Axis -- Inclination


ASolar=0.58;                                                            % SunLit Surface Area in m^2.
QSolar=0.05;
Cps=0.2866;                                        % Center of Solar Pressure in meters.
Cpa=0.658;                                         % Center of AeroDynamic Pressure in meters. - Max. Length of Deck=1.974meters.
Cms=0;                                             % Center of Mass in meters.    
AIS=degtorad(0);                                   % Angle of Incidence in radians
SDM=[-1;1;1];                                      % Spacecraft residual Dipole Moment A.m^2
Mag=7.943*10^15;                                   % Earth Magnetic Constant in Tesla.
DragC=2;                                           % Drag Coeffiecient.
RArea=1;                                           % Ramp Area in m^2.
%Adens in Kg/m^3   Converted from g/cm-3  by Multiplying by 10^3;
Adens700=2.503e-14; Adens750=1.381e-14; Adens800=8.353e-15;  Adens850=5.534e-15;  Adens900=3.962e-15;  Adens950=3.006e-15;  Adens1000=2.373e-15;
%-------------------------------------------------------------------------%
CurrentTA = OrbitElements(6)*pi/180 + TAdot*t;  % Radians 
RadiusPerigee=(SMA*(1-ecc^2))/(1+ecc*cos(CurrentTA));                             % Norm of Radius From the Center of Earth to S/C in Orbit.
OrbitVelocity=sqrt(GConstant*((2/RadiusPerigee)-(1/SMA)))*10^3;                 % Norm of Velocity in m/sec.
Alpha3=-(AOP+CurrentTA);                                                          % TAdot*t= Gives True Anomaly. Alpha3 is AOP + TrueAnomaly.
RotAOP=[cos(Alpha3) sin(Alpha3) 0;-sin(Alpha3) cos(Alpha3) 0;0 0 1];            % Z-Axis Rotation Around True Anomaly.
RUV=RotRANN*RotINC*RotAOP*[1;0;0];                                              % Unit Vector in ECI Frame. - RTN TO ECI Frame.
RVectorECI=RotRANN*RotINC*RotAOP*[RadiusPerigee;0;0];                           % Radius Vector in ECI Frame.
%Quaternion Transformation
QMatrix=[qb(1)^2-qb(2)^2-qb(3)^2+qb(4)^2 2*(qb(1)*qb(2)+qb(3)*qb(4)) 2*(qb(1)*qb(3)-qb(2)*qb(4));2*(qb(1)*qb(2)-qb(3)*qb(4)) -qb(1)^2+qb(2)^2-qb(3)^2+qb(4)^2 2*(qb(2)*qb(3)+qb(1)*qb(4));2*(qb(1)*qb(3)+qb(2)*qb(4)) 2*(qb(2)*qb(3)-qb(1)*qb(4)) -qb(1)^2-qb(2)^2+qb(3)^2+qb(4)^2];
%--------------------ECI To ECEF To Lat Long-------------------%
if(strcmp(Planet,'Eart'))   % Magentic Field Only For Earth
% Julian Date & Rvector in ECI -O/P R Vector in ECEF
THETAEarth = JD2GAST(JD);              % ECI to ECEF - Approx. Rotation around Z-Axis of ECI to go to ECEF. 
RVectorECEF=[cosd(THETAEarth) sind(THETAEarth) 0;-sind(THETAEarth) cosd(THETAEarth) 0;0 0 1]*[RVectorECI(1);RVectorECI(2);RVectorECI(3)];                 

% Rotation from ECEF to Magnetic Frame, 11.5 Around Z-Axis, -19.3 Around Y-Axis, 
RzECEF2Mag=[cosd(11.5) sind(11.5) 0;-sind(11.5) cosd(11.5) 0;0 0 1];
RyECEF2Mag=[cosd(-19.3) 0 -sind(-19.3);0 1 0;sind(-19.3) 0 cosd(-19.3)];  

% SpaceCraft Position in Magnetic Frame
RVectorMagnetic=RyECEF2Mag*RzECEF2Mag*RVectorECEF;

[LongMag,LatMag]=cart2sph(RVectorMagnetic(1),RVectorMagnetic(2),RVectorMagnetic(3));
LongMag=LongMag*180/pi;
LatMag=LatMag*180/pi;

%Magnetic Field in Magnetic Frame
BMagnetic=-(Mag/(RadiusPerigee*10^3)^3)*[3*sin(LatMag)*cos(LatMag)*cos(LongMag);3*sin(LatMag)*cos(LatMag)*sin(LongMag);3*sin(LatMag)^2-1];                                                        


SMDECI=QMatrix'*SDM;        % Body To Inertial, SpaceCraft Residual Dipole. Amp-m^2.
SMDECEF=[cosd(THETAEarth) sind(THETAEarth) 0;-sind(THETAEarth) cosd(THETAEarth) 0;0 0 1]*SMDECI;
SMDMag=RyECEF2Mag*RzECEF2Mag*SMDECEF;

% Magnetic Torque
TMagneticx=SMDMag(2)*BMagnetic(3)-SMDMag(3)*BMagnetic(2);
TMagneticy=-SMDMag(1)*BMagnetic(3)+SMDMag(3)*BMagnetic(1);              
TMagneticz=SMDMag(1)*BMagnetic(2)-SMDMag(2)*BMagnetic(1);
end
dca=QMatrix*RUV;           % Inertial to Body. Radius Vector in Body Frame. 
dcx=dca(1); dcy=dca(2); dcz=dca(3);

%Gravity Gradient Torque
TGx=((3*GConstant)/RadiusPerigee^3)*((Isc(3,3)-Isc(2,2))*dcy*dcz+Isc(2,3)*(dcy^2-dcz^2)+Isc(1,3)*dcx*dcy-Isc(1,2)*dcx*dcz);
TGy=((3*GConstant)/RadiusPerigee^3)*((Isc(1,1)-Isc(3,3))*dcz*dcx+Isc(3,1)*(dcz^2-dcx^2)+Isc(2,1)*dcy*dcz-Isc(2,3)*dcy*dcx);
TGz=((3*GConstant)/RadiusPerigee^3)*((Isc(2,2)-Isc(1,1))*dcx*dcy+Isc(1,2)*(dcx^2-dcy^2)+Isc(3,2)*dcz*dcx-Isc(3,1)*dcz*dcy);

% Solar Radiation Torque
% Solar Panel Sizes Middle Panel - a=662mm, b=880mm, h=889mm.
% Area=0.6854m^2, Solar Panel CG from Model Frame=[880,-0.31,662]mm
% Normal Vector in Body Frame = BPanel, Az=180,El=15
% Solar Panel Sizes Side1 on - Y Axis Panels - a=428mm, b=735mm, h=890mm.
% Area=0.5175m^2, Solar Panel CG from Model Frame=[671,-599,646]mm
% Normal Vector in Body Frame = APanel, Az=130,El=15
% Solar Panel Sizes Side1 on + Y Axis Panels - a=428mm, b=735mm, h=890mm.
% Area=0.5175m^2, Solar Panel CG from Model Frame=[671,599,646]mm
% Normal Vector in Body Frame = CPanel, Az=230,El=15

% As the Solar Panel is on Negative X-Axis of Spacecraft. 
% Torque Contribution from Middle Solar Panel
AMiddle=0.6854; % m^2.
ASide1=0.5175;
ASide2=0.5175;
SunVectorBody=QMatrix*SunVectorJ2000;   % Solar Panel on -X-Axis. This Assumption Assumes when qb=[0;0;0;1]. The Spacecraft is Pointing towards the Sun. 
AngleMiddlePanel=acosd(dot(SunVectorBody,BPanel));
AngleSide1Panel=acosd(dot(SunVectorBody,APanel));
AngleSide2Panel=acosd(dot(SunVectorBody,CPanel));
RMiddle = [0.88;-0.00031;0.662] - rcg;
RSide1 = [0.671;-0.599;0.646] - rcg;
RSide2 = [0.671;0.599;0.646] - rcg;
    TPSolarForce = (1366/(3*10^8))*0.21;  
if(AngleMiddlePanel<90)
    SolarForceMiddle = TPSolarForce*cosd(AngleMiddlePanel)*AMiddle;
    ForceMiddleComponents = SunVectorJ2000*SolarForceMiddle;
    TMiddle = cross(RMiddle,ForceMiddleComponents);
else
    TMiddle=zeros(3,1);
end
if(AngleSide1Panel<90)
    SolarForceSide1 = TPSolarForce*cosd(AngleSide2Panel)*ASide1;
    ForceSide1Components = SunVectorJ2000*SolarForceSide1;
    TSide1=cross(RSide1,ForceSide1Components);
else
    TSide1=zeros(3,1);
end
if(AngleSide2Panel<90)
    SolarForceSide2 = TPSolarForce*cosd(AngleSide1Panel)*ASide2;
    ForceSide2Components = SunVectorJ2000*SolarForceSide2;
    TSide2=cross(RSide2,ForceSide2Components);
else
    TSide2=zeros(3,1);
end
Tsolarx = TMiddle(1) + TSide1(1) + TSide2(1);
Tsolary = TMiddle(2) + TSide1(2) + TSide2(2);
Tsolarz = TMiddle(3) + TSide1(3) + TSide2(3);
% Tsolarx=(1366/(3*10^8))*AMiddle*(1+QSolar)*(Cps-Cms)*cos(AIS);

Height=RadiusPerigee-PlanetRadius;
if(Height>0 && Height<700)
    Adens=Adens700;
elseif(Height>700 && Height<750)
    Adens=Adens750;
elseif(Height>750 && Height<800)
    Adens=Adens800;
elseif(Height>800 && Height<850)
    Adens=Adens850;
elseif(Height>850 && Height<900)
    Adens=Adens900;
elseif(Height>900 && Height<950)
     Adens=Adens950;
elseif(Height>950 && Height<1000)
     Adens=Adens1000;
elseif(Height>=1000)
    Adens=0;
end

% Atmospheric Drag
TAtmos=0.5*Adens*DragC*RArea*OrbitVelocity^2*(Cpa-Cms);                 
Text=zeros(3,1);
if(strcmp(Planet,'Moon'))
    SolarTorque = [Tsolarx;Tsolary;Tsolarz];
    MagTorque = zeros(3,1);
    GravTorque= [TGx;TGy;TGz];
    Text(1)=TGx+Tsolarx;
    Text(2)=TGy+Tsolary;
    Text(3)=TGz+Tsolarz;
elseif(strcmp(Planet,'Eart'))
    SolarTorque = [Tsolarx;Tsolary;Tsolarz];
    MagTorque = [TMagneticx;TMagneticy;TMagneticz];
    GravTorque= [TGx;TGy;TGz];
    Text(1)=TGx+Tsolarx+TMagneticx+TAtmos;
    Text(2)=TGy+TMagneticy+Tsolary;
    Text(3)=TGz+TMagneticz+Tsolarz;
end    
end    
