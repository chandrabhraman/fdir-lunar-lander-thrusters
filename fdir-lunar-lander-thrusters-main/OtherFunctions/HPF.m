function [Output]=HPF(Input,PrevInput,PrevOutput,tinc,f2)
% Ver1. Created 29-7-2015

% High Pass Filter
%%% Input = Current Input to Filter, Previous Input to Filter, Previous
%%% Output of Filter, StepSize, F2. 
% For Derivation Refer. PR-Modulator Document
Output=(Input-PrevInput)/(1+f2*tinc)+PrevOutput/(1+f2*tinc);

end