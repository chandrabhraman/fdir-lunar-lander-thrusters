
%% This part of the code plots graphs for (measured acceleration - Ti as fault mode acceleration' )and stores in folder '
%{
plot(T(1).debug(:,1), T(1).debug(:,2))
hold on
plot(T(2).debug(:,1), T(2).debug(:,2))
hold on
plot(T(3).debug(:,1), T(3).debug(:,2));
hold on;
plot(T(4).debug(:,1), T(4).debug(:,2));
legend('T1','T4','T5','T8')
xlabel('time (seconds)')
ylabel('Angular acc (deg/s^2)')
title('plot for alpha_x with T5 as failure')

plot(T(1).debug(:,1), T(1).debug(:,3))
hold on
plot(T(2).debug(:,1), T(2).debug(:,3))
hold on
plot(T(3).debug(:,1), T(3).debug(:,3));
hold on;
plot(T(4).debug(:,1), T(4).debug(:,3));
legend('T1','T4','T5','T8')
xlabel('time (seconds)')
ylabel('Angular acc (deg/s^2)')
title('plot for alpha_y with T5 as failure')

plot(T(1).debug(:,1), T(1).debug(:,4))
hold on
plot(T(2).debug(:,1), T(2).debug(:,4))
hold on
plot(T(3).debug(:,1), T(3).debug(:,4));
hold on;
plot(T(4).debug(:,1), T(4).debug(:,4));
legend('T1','T4','T5','T8')
xlabel('time (seconds)')
ylabel('Angular acc (deg/s^2)')
title('plot for alpha_z with T5 as failure')

plot(T(1).debug(:,1), T(1).debug(:,5))
hold on
plot(T(2).debug(:,1), T(2).debug(:,5))
hold on
plot(T(3).debug(:,1), T(3).debug(:,5));
hold on;
plot(T(4).debug(:,1), T(4).debug(:,5));
legend('T1','T4','T5','T8')
xlabel('time (seconds)')
ylabel('Linear acc (m/s^2)')
title('plot for a_x with T5 as failure')

plot(T(1).debug(:,1), T(1).debug(:,6))
hold on
plot(T(2).debug(:,1), T(2).debug(:,6))
hold on
plot(T(3).debug(:,1), T(3).debug(:,6));
hold on;
plot(T(4).debug(:,1), T(4).debug(:,6));
legend('T1','T4','T5','T8')
xlabel('time (seconds)')
ylabel('Linear acc (m/s^2)')
title('plot for a_y with T5 as failure')

plot(T(1).debug(:,1), T(1).debug(:,7))
hold on
plot(T(2).debug(:,1), T(2).debug(:,7))
hold on
plot(T(3).debug(:,1), T(3).debug(:,7));
hold on;
plot(T(4).debug(:,1), T(4).debug(:,7));
legend('T1','T4','T5','T8')
xlabel('time (seconds)')
ylabel('Linear acc (m/s^2)')
title('plot for a_z with T5 as failure')

%}

%% This part of the code collects telemetry data & then runs Thruster_FDIR'
%% Second part plots graphs and stores in folder ' To store the figures of all failures modes.
close all
clc
clear all
SSU_Gyro_FDIR
close all
Thruster_FDIR
%{
cd './LinearAcceleration/Failure plots T8/Graphs'
saveas(figure(1),'FiringHistory.jpg');
saveas(figure(2),'ActiveAngularAcc.jpg');
saveas(figure(3),'ActiveLinearAcc.jpg');
saveas(figure(4),'InactiveAngularAcc.jpg');
saveas(figure(5),'InactiveLinearAcc.jpg');
save('ThrusterFDIR.mat','ThrusterFDIR')
save('T.mat','T')
cd '../../../'
%}
