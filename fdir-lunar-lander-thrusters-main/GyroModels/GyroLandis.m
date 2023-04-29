% Gyro Model from Landis Markley, Pg:-147

function [WMeas,Bias]=GyroLandis(WInput,BiasPrev,StepSize,RandomNumber)
%%%
% Input = Body Rates without Noise, CurrentBias, StepSize
% Output = Corrupted Rates, UpdatedBias.
%%%
%%%
% Date, Changes, Made by Name
%Ver1 - 15/8/15 - Random Number Function Kept in Main
%%%
ARW=6.399540590645874e-07;     %rad/sec^1/2 Angular Random Walk
RRW=2.278624301214819e-10;     %rad/sec^3/2 Rate Random Walk

BiasNoise=RRW*sqrt(StepSize)*RandomNumber(1);
Bias=BiasPrev+BiasNoise;
Noise=sqrt(ARW^2/StepSize+(1/12)*RRW^2*StepSize)*RandomNumber(2);
WMeas=WInput+0.5*(Bias+BiasPrev)+Noise;

end
