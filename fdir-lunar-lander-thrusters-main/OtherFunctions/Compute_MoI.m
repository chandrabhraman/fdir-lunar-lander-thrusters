function [MoI_total_CGFrame] = Compute_MoI(r_cg,theta_f,theta_ox,m_f,m_ox,RadiusTank,Fuel_Dens,Oxi_Dens,Fuel_Tank_Position,Oxi_Tank_Position,DryMass)


% MoI_dry_StructureFrame=[91.66 -2.62 1.34;
%                         -2.62 61.62 3.10;
%                         1.34 3.10 94.02];
% 
% 
% Stru2GNC=[1 0 0;
%           0 0 -1;
%           0 1 0];
% 
% MoI_dry_ModelFrame= Stru2GNC*MoI_dry_StructureFrame*Stru2GNC';

MoI_dry_ModelFrame=[136.107196 -1.478583  -1.249036;1.478583   128.413548  11.916741;-1.249036   11.916741    93.642413];

% Compute the Integrals Associated with Fuel and Oxidizer Levels
[C1_ox,C2_ox,C3_ox,C4_ox] = Compute_liquid_level_integrals(theta_ox);       
[C1_f,C2_f,C3_f,C4_f] = Compute_liquid_level_integrals(theta_f);

% MOI Matrix contribution from oxidizer, given in model frame
% Get coordinates of center of oxidizer tank, in model frame
x0_ox = Oxi_Tank_Position(1); 
y0_ox = Oxi_Tank_Position(2); 
z0_ox = Oxi_Tank_Position(3); 
MoI_ox_ModelFrame = Compute_MoI_terms(Oxi_Dens,RadiusTank,C1_ox,C2_ox,C3_ox,C4_ox,x0_ox,y0_ox,z0_ox);      

% MoI matrix contribution from fuel, given in model frame
% Get coordinates of center of fuel tank, in model frame
x0_f = Fuel_Tank_Position(1);
y0_f = Fuel_Tank_Position(2);
z0_f = Fuel_Tank_Position(3);
MoI_f_ModelFrame = Compute_MoI_terms(Fuel_Dens,RadiusTank,C1_f,C2_f,C3_f,C4_f,x0_f,y0_f,z0_f);

% Compute Total MOI matrix about the model frame
MoI_total_ModelFrame = MoI_dry_ModelFrame + MoI_ox_ModelFrame + MoI_f_ModelFrame;

% Compute total MoI matrix about (instantaneous) CG frame (We assume that
% the frame at CG is parallel to the model frame). We use the formula from
% Likins, p422.

% Compute total mass of spacecraft at current time
m_total = DryMass + m_ox + m_f;              


% Note : 'r_cg' is taken to be From the model frame origin to the current cg
%location, expressed in  the model frame.
% SkewSym
M = [0 -r_cg(3) r_cg(2);r_cg(3) 0 -r_cg(1);-r_cg(2) r_cg(1) 0];
    
MoI_total_CGFrame = MoI_total_ModelFrame + m_total*M*M;

end