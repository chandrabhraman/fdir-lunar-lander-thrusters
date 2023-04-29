function [PosVelDot]=FunPosVel(t,PosVelo,Inertial_Acc)
global JulianArray SineCoef CosCoef Re GM DeltaAT DeltaUTC PolarX PolarY Degree DeltaPsi DeltaEpsilon

JDT=JulianArray+t/86400; %Convert present time to julian data with increment from JulianArray
Rx=PosVelo(1);
Ry=PosVelo(2);
Rz=PosVelo(3);
Vx=PosVelo(4);
Vy=PosVelo(5);
Vz=PosVelo(6);

PosInertial=[Rx;Ry;Rz];

et=cspice_str2et(sprintf('JD%0.9f',JDT)); %Get ephemeris Time from SPICE function
DCMJ2000toFixed=FK5(JDT,DeltaAT,DeltaUTC,[DeltaPsi DeltaEpsilon],[PolarX PolarY]); % ECI2ECEF
[psun]=cspice_spkpos('sun',et,'J2000','NONE','earth');
[pmoon]=cspice_spkpos('moon',et,'J2000','NONE','earth');
PointMassAcc=PointMassGravityEarth(PosInertial,psun,pmoon);     % Gravity Due to Moon and Sun

NSG=NonSphericalGravityV1(PosInertial,DCMJ2000toFixed,SineCoef,CosCoef,Degree,Re,GM);   % NonSpherical Gravity
TotalAcc = NSG+PointMassAcc+Inertial_Acc;

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
