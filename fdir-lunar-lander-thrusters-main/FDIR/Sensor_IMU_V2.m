% Three axis Accelerometer model - [30-04-2016]

function Ameas=Sensor_IMU_V2(wb, wb_dot, d, IMU_Accelerometer_SF_MA, Acc_true, tinc)


% d = [(Accelerometer_Location(1) -  r_cg(1));
%      (Accelerometer_Location(2) -  r_cg(2));
%      (Accelerometer_Location(3) -  r_cg(3))];

Aimeas = Acc_true + cross(wb,cross(wb,d)) + cross(wb_dot,d);

accelero_sig_VRW=(0.0025)*(1/3.28)*(1/60);      % Velocity random walk (m/s) per rt(sec)
accelero_sig_beta=accelero_sig_VRW*sqrt(tinc);  % !! IMPORTANT !! The integration time step is the same as the sensor sampling time
accelero_beta=accelero_sig_beta*randn(3,1);     % The beta measurement error signal generated

% Accelerometer Model bias signal
accelero_sig_bias=5*1e-6*9.81;                  % 5 micro g  (m/s2) 1 sigma
accelero_bias=accelero_sig_bias*randn(3,1);     % Moving Bias signal generated (m/s)

Ameas = (eye(3,3)+ IMU_Accelerometer_SF_MA)*Aimeas + accelero_bias + accelero_beta;

end
