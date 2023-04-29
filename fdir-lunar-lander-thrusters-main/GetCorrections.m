function [DeltaAT,DeltaUTC,PolarX,PolarY,DeltaPsi,DeltaEpsilon]=GetCorrections(JD)

load taiutc.mat
load eopc04.mat
MJD = JD - 2400000.5;   % Modified Julian Date for Getting Polar Motion and DeltaUTC
[TaiRow,~]=size(taiutc);    % Get taiutc Row Size
[eopRow,~]=size(eopc04);

% Trying to find the Closest Match for TAI-UTC , UTC1-UTC, Polar Motion,
% dPsi and dEpsilon - Not Interpolating . - Data Set Taken From GMAT -
% Which is Coarse- Need to replace with a finer one. - 26/10/15
for i = 1:TaiRow
    Store1(i) = (taiutc(i,1) - JD)^2;
end

for i=1:eopRow
    Store2(i) = (eopc04(i,4) - MJD)^2;
end

[~,TaiIndex]=min(Store1);
[~,EopIndex]=min(Store2);

DeltaAT = taiutc(TaiIndex,2) + (MJD - taiutc(TaiIndex,3))*taiutc(TaiIndex,4);
DeltaUTC = eopc04(EopIndex,7);
PolarX = (eopc04(EopIndex,5)/3600)*pi/180;
PolarY = (eopc04(EopIndex,6)/3600)*pi/180;
DeltaPsi = (eopc04(EopIndex,9)/3600)*pi/180;
DeltaEpsilon = (eopc04(EopIndex,10)/3600)*pi/180;
end