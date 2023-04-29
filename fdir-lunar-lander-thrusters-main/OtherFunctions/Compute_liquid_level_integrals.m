function [C1,C2,C3,C4] = Compute_liquid_level_integrals(theta)

% This function computes the integrals associated with oxidizer/fuel
% levels.

C1 = (1/80)*cos(5*theta) - (5/48)*cos(3*theta) + (5/8)*cos(theta) + 8/15;
C2 = (3/4)*cos(theta) - (1/12)*cos(3*theta) + 2/3;
C3 = (1/48)*cos(3*theta) - (1/80)*cos(5*theta) + (1/8)*cos(theta) + 2/15;
C4 = -(1/4)*(sin(theta))^4;

end