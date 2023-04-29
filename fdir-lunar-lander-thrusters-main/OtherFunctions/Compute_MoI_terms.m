function [MoI_ModelFrame] = Compute_MoI_terms(rho,R,C1,C2,C3,C4,x0,y0,z0)

% For either Fuel/Oxidizer, this function computes the MoI matrix about the
% model frame.

Ixx_ModelFrame = (0.25*rho*pi*R^5)*C1 + (rho*pi*R^3*y0^2)*C2 ...
                    + (rho*pi*R^3*z0^2)*C2 + (rho*pi*R^5)*C3 ...
                    + (2*rho*pi*R^4*z0)*C4;

Iyy_ModelFrame = (0.25*rho*pi*R^5)*C1 + (rho*pi*R^3*x0^2)*C2 ...
                    + (rho*pi*R^3*z0^2)*C2 + (rho*pi*R^5)*C3 ...
                    + (2*rho*pi*R^4*z0)*C4;
                
Izz_ModelFrame = (0.5*rho*pi*R^5)*C1 + (rho*pi*R^3*y0^2)*C2 ...
                    + (rho*pi*R^3*x0^2)*C2;

Ixy_ModelFrame = -(rho*pi*R^3*x0*y0)*C2;
Iyx_ModelFrame = Ixy_ModelFrame;

Iyz_ModelFrame = -(rho*pi*R^3*y0*z0)*C2 - (rho*pi*R^4*y0)*C4;
Izy_ModelFrame = Iyz_ModelFrame;

Izx_ModelFrame = -(rho*pi*R^3*z0*x0)*C2 - (rho*pi*R^4*x0)*C4;
Ixz_ModelFrame = Izx_ModelFrame;

MoI_ModelFrame = [
                Ixx_ModelFrame Ixy_ModelFrame Ixz_ModelFrame;
                Iyx_ModelFrame Iyy_ModelFrame Iyz_ModelFrame;
                Izx_ModelFrame Izy_ModelFrame Izz_ModelFrame
               ];
end