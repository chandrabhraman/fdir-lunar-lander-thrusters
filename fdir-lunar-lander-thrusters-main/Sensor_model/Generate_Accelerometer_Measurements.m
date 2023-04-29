function Ameasured = Generate_Accelerometer_Measurements(wb, acc_true, wb_dot)

% 01 Jun 2016

close all
clc

N = size(wb,1); % wb is true body rates, [Nx3] -> should be equal to the number of true acceleration samples
h = 0.016; % Simulation step size
Accelerometer_Location = [-772.8; -232.9032; 321.5598] *1e-3; % [m]
r_cg = [-3.3; 0; 555.9] *1e-3; %[m]

% Generate bias signal
IMU_Accelerometer_bias = 5*1e-6*9.81 * randn(N,3); % [Nx3]

% Generate beta signal
IMU_Accelerometer_beta = (0.0025*1/3.28*1/60) * sqrt(h) *  randn(N,3); % [Nx3]

% IMU SF and MA matrix
IMU_Accelerometer_SF_MA = [ 0.0001                  8.29e-05                  8.29e-05
                            8.29e-05                    0.0001                  8.29e-05
                            8.29e-05                  8.29e-05                    0.0001];

for i = 1 : N
    current_wb = wb(i,:)';  % [3x1]
    current_acc_true = acc_true(i,:)'; % [3x1]
    current_wb_dot = wb_dot(i,:)'; % [3x1]
    Ameasured(i,:) = Sensor_IMU_V2(i, current_wb, Accelerometer_Location, r_cg, ...
                                   IMU_Accelerometer_SF_MA, IMU_Accelerometer_bias, ...
                                   IMU_Accelerometer_beta, current_acc_true, current_wb_dot);
end


end