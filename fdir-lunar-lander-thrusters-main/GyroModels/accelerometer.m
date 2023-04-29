
function nga_meas=accelerometer(nga,tinc)
% Use Accelerometer Model  
% Input - NGA= m/sec^2

accelero_sig_VRW=(0.0025)*(1/3.28)*(1/60);      % Velocity random walk (m/s) per rt(sec)
accelero_sig_beta=accelero_sig_VRW*sqrt(tinc);  % !! IMPORTANT !! The integration time step is the same as the sensor sampling time
accelero_beta=accelero_sig_beta*randn(3,1);     % The beta measurement error signal generated

% Accelerometer Model bias signal
accelero_sig_bias=5*1e-6*9.81;                  % 5 micro g  (m/s2) 1 sigma
accelero_bias=accelero_sig_bias*randn(3,1);     % Moving Bias signal generated (m/s)
nga_meas=nga+accelero_bias+accelero_beta; 

end

