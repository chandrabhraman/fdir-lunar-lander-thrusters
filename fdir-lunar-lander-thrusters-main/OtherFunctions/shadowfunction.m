% Shadow Function - From GMAT Math Spec. 

function [SRP,Eclipse]=shadowfunction(SatPos,CB2Sun,CBRadius,mass_sc)

% SatPos - Position of Statellite from Central Body  - Km - 3x1 Vector
% CB2Sun - Position of Sun from Central Body    - Km - 3x1 Vector
% mass_sc - Mass of Spacecraft - KG

if(isrow(SatPos))
    SatPos=SatPos';
elseif(isrow(CB2Sun))
    CB2Sun=CB2Sun';
end

RS = 0; % Sun Radius
RB = CBRadius;  % Central Body radius - Km
SC2Sun = CB2Sun - SatPos;   % Position of S/C to Sun - Km

RSunDash = asin(RS/norm(SC2Sun));   % Apparent Radius of Sun
RCBDash = asin(RB/norm(SatPos));    % Apparent Radius of CB

DDash = acos(-SatPos'*SC2Sun/(norm(SatPos)*norm(SC2Sun)));

if(DDash>=(RSunDash+RCBDash))
    Eclipse=0;  % No Eclipse 
elseif(DDash<=(RCBDash-RSunDash))
    Eclipse=100;    % Full Shadow - Eclipse 
elseif(DDash<(RSunDash+RCBDash) && DDash>abs(RSunDash-RCBDash))  %!! Verify Equation Term Mismatch RSunDash - PDFNo: 69 - GMAT Math Spec.
    c1= (DDash^2+RSunDash^2-RCBDash^2)/(2*DDash);
    c2=sqrt(RSunDash^2-c1^2);
    A= RSunDash^2*acos(c1/RSunDash) + RCBDash^2*acos((DDash-c1)/RCBDash) - DDash*c2;
    Eclipse =  100*A/(pi*RSunDash^2);
else
    Eclipse = 100*(RCBDash^2/RSunDash^2);
end

UnitSC2Sun = SC2Sun/norm(SC2Sun);

PercSun=(100-Eclipse)/100;
SolarFlux = 1367; % W/m^2
Absorption = 0.21;  % !! Need to Find out
Cr = 1 + Absorption;
SL = 3 * 10^8;  % m/Sec
AU = 149597870700/1000; % In Km's
RAu = norm(SC2Sun)/AU;  % In AU   
Area = 1;   % m^2 !! Need to Find Out.
SRP = -Cr*(Area/mass_sc)*(PercSun*SolarFlux/SL)*(1/RAu^2)*UnitSC2Sun;

% Cr - Coefficient of Reflectivity 
% Area - Aera of Spacecraft 
% PercSun - Percentage of Sun Visbile to SpaceCraft
% SolarFlux - SolarFlux at 1 AU
% SL - Speed of Light 
% RAu - Distance from S/C to Sun in AU.
% UnitSC2Sun - Position Vector of Sun as Seen from Spacecraft. 

end
    
    
    

