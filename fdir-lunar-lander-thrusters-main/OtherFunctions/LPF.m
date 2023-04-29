function [Output]=LPF(Input,PrevOutput,tinc,f2)
% Ver1. Created - 29/7/2015

% Low pass Filter
%%% Input - Current Input, Previous Filter Output, Step Size, f1, f2
% For Derivaiton Refer PR Modulator Document
Output=(tinc/(1+f2*tinc))*Input+PrevOutput*(1/(1+f2*tinc));

end