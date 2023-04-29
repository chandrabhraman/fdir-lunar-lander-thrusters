function y=RK4(h,t,x,Inertial_Acc)

K1 = FunPosVel(t,x,Inertial_Acc);
K2 = FunPosVel(t+h/2,x+(h/2)*K1,Inertial_Acc);
K3 = FunPosVel(t+h/2,x+(h/2)*K2,Inertial_Acc);
K4 = FunPosVel(t+h,x+h*K3,Inertial_Acc);

y = x+(h/6)*(K1+2*K2+2*K3+K4);
end
