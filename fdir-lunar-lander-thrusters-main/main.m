% Normal Mode
clear all
close all
clc

%% Setup
addpath(genpath('Mice_data'));
addpath(genpath('OtherFunctions'));
addpath(genpath('FunctionsKinematicDynamics'));
addpath(genpath('DisturbanceTorqueModel'));
addpath(genpath('GyroModels'));

% Load Mice Kernels
cspice_furnsh('Mice_data\MICE\kernels\lsk\naif0011.tls.pc')
cspice_furnsh('Mice_data\MICE\kernels\pck\pck00010.tpc')
cspice_furnsh('Mice_data\MICE\kernels\pck\earth_fixed.tf')
cspice_furnsh('Mice_data\MICE\kernels\pck\earth_720101_070426.bpc')
cspice_furnsh('Mice_data\MICE\kernels\spk\de421Planets.bsp')
cspice_furnsh('Mice_data\MICE\kernels\pck\moon_pa_de421_1900-2050.bpc')
cspice_furnsh('Mice_data\MICE\kernels\pck\moon_080317.tf')

load OrbitElements.mat
load Block1.mat
load ThrusterParameters.mat
load PWPFParameters.mat
load PositionCanting.mat
load FuelTankParameters.mat
load C2_tank.mat

%% Simulation Parameters

tinc=0.016; % Time Step in Seconds
tfF=round(30/tinc);
tinc16 = 0.016;

%% Environment parameters
global JulianArray SineCoef CosCoef Re GM DeltaAT DeltaUTC PolarX PolarY Degree DeltaPsi DeltaEpsilon NutationData
JulianArray = 2457956.2948392;
CurrentOrbit = 'G1';
OrbitElements = [eval(sprintf(strcat('Orbit','.',CurrentOrbit,'.','Ecc')));eval(sprintf(strcat('Orbit','.',CurrentOrbit,'.','SMA')));eval(sprintf(strcat('Orbit','.',CurrentOrbit,'.','Inc')));eval(sprintf(strcat('Orbit','.',CurrentOrbit,'.','RAAN')));eval(sprintf(strcat('Orbit','.',CurrentOrbit,'.','AOP')))]';
RefAttitude=eval(sprintf(strcat('Orbit','.',CurrentOrbit,'.','Attitude')));
[DeltaAT,DeltaUTC,PolarX,PolarY,DeltaPsi,DeltaEpsilon]=GetCorrections(JulianArray);
load NutationData.mat
load EGM96.mat
Degree = 30;
g0=GravityEarth;

%% Normal-Mode Parameters

% First Burn Position and Velocity Vector data (data is @ every 16ms)
% [Cols 1-3 are ECI positions in Km, Cols 4-6 are ECI velocities in Km/s]
% PosVel = csvread('ENB Position and velocity.csv',4,1);
mu = 398600;
oev = [OrbitElements(2);OrbitElements(1);OrbitElements(3);OrbitElements(4);OrbitElements(5);0];
[Pos,Vel] = orb2eci(mu, oev);
PosVel = [Pos;Vel]';

% Gyro Bias
Bias=(1*pi)/(180*3600);     % rad/sec
Bias=[Bias;Bias;Bias];
tau=0.1;                    % time Constant for the low pass filter
w_lpf=[0.0;0.0;0.0];

%-----------------------Filter-------------------------------------------%
TKG=[0.0442998372043458,-5.19682084201292e-08,-2.10170349459057e-06,-5.19682084441953e-08,0.0442998306257591,1.96283505768552e-06,-0.000134509023653650,0.000125621443691780,0.00688656615407449,-0.000439166878611658,1.93661058018496e-07,-2.44250350854836e-07,-1.97857408875370e-07,-0.000439167317449383,2.30932156991861e-07,1.55320813446343e-06,-1.45107362369680e-06,-6.82584906927437e-05];
KG=[TKG(1:3);TKG(4:6);TKG(7:9);TKG(10:12);TKG(13:15);TKG(16:18)];
Starnoise=[1/3600 1/3600 8/3600];

% Estimates
QEstimate=qnorm(RefAttitude');
wbEstimate=[0.1*pi/180;0.1*pi/180;0.1*pi/180]; % Angular velocity in bf
BiasEstimate=[0;0;0];
RandomNumbers=randn(5,tfF);

%% Spacecraft Parameters for simulation
Initial_Fuel_mass           = 212.9;  % Kg
Initial_Oxi_mass            = 187.23; % Kg
DryMass                     = 199.87; % Kg
Dry_CG                      = [10; 0; 443.61]/1000; % m
Initial_SC_mass             = Initial_Fuel_mass + Initial_Oxi_mass + DryMass; % Kg
Main_Engine_ThrustElevation = 89.8*pi/180;  % rad
Main_Engine_ThrustAzimuth   = 0*pi/180;     % rad
Dry_MoI_BodyFrame = [ 136.107196 -1.478583  -1.249036;
    -1.478583   128.413548  11.916741;
    -1.249036   11.916741    93.642413]; % MoI of dry spacecraft wrt Body Frame, Kg-m^2

RCG_guess = [ -0.00334868333333335
    0
    0.555904783422407]; % Center of mass of spacecraft @ Launch

IMsc = [ 120.558405695916     -1.478583          8.78189256276734
    -1.478583          191.753164043163   11.916741
    8.78189256276734    11.916741        196.19658148162]; % Kg-m^2

InvIMsc = inv(IMsc);

NominalThrust = ACTThrustSunPointing; % N (Nominal thrust of attitude control thrusters in ON modulation)
ISP           = ACTISPSunPointing;    % s
m_dot         = NominalThrust/(ISP*g0); % Kg/s
m_dot_LAM     = LAMThrust/(LAMISP*g0);  % Kg/s
LAM_Status = 0; % Switched OFF during Normal Mode

% Main Engine Unit Thrust Vector in the spacecraft bodyframe
[CartX,CartY,CartZ]=sph2cart(Main_Engine_ThrustAzimuth, Main_Engine_ThrustElevation, 1);

% ------                Solar Panel Placements --------------------------%
[BPanel(1),BPanel(2),BPanel(3)]=sph2cart(180*pi/180,15*pi/180,1);  % Panel B
[APanel(1),APanel(2),APanel(3)]=sph2cart(130*pi/180,15*pi/180,1);  % Panel A
[CPanel(1),CPanel(2),CPanel(3)]=sph2cart(230*pi/180,15*pi/180,1);  % Panel C

%% Control and PWPFM Parameters
Kp = [8.3815176359162        -0.102492803529698        0.731825607944399
    0.0690572743520295      8.92913451543557         0.616298781133795
    0.110564383584068       0.139027756972378        2.06085519165111];

Kd = [35.905368325061        -0.439066291005856           3.1350489427304
    0.381243402312288          49.2949316970645          3.40239093368531
    0.885844268352307          1.11389344075903          16.5116170347943];

Uon  = PWPF16Uon;
Uoff = PWPF16Uoff;
Tm   = PWPF16Tm;
Km   = PWPF16Km;

wref =[0.0*pi/180;0.0*pi/180;0.0*pi/180]; % Ref body rate(rad/s)

TLimit = [11.2800041354447          15.4004236416049          2.70891157160411]; % Torque limits used to normalize and saturate commanded torques
qeLimit =[1.34581881533101 1.72473867595819 1.31445993031359]';
%% Initial conditions

% Initial time
t0 = 0;

% Initial Position and Velocity
r0 = PosVel(1,1:3)';     % Position Vector in Inertial [Km]
v0 = PosVel(1,4:6)';     % Velocity Vector in Inertial [Km/s]

% Initial Attitude
q0 = qnorm(RefAttitude'); % Normalizing Quaternions

% Initial body rates
w0=[0.1*pi/180;0.1*pi/180;0.1*pi/180];

%% Initialize stuff
Fuel_Mass = Initial_Fuel_mass;
Oxi_Mass = Initial_Oxi_mass;
DeltaV   = 0;
Text     = zeros(3,1);
mass_spent_LAM = 0;
mass_spent_ACT = 0;
mass_spent     = 0;
ThrustFilterPrevious = zeros(TotalThrusters,1);
py = [0 0 0]';
Um = [0 0 0]';
LP = [0 0 0]';
LNP = [0 0 0]';

%% Storage
Output.States   = zeros(tfF+1, 14);
Output.SC_Param = zeros(tfF+1, 7);
Output.ACT_CMD = zeros(tfF+1, 5);
Output.Control = zeros(tfF+1,5);
Output.Estimate = zeros(tfF+1,10);

%% Loop
for t=1:tfF

    %---------------------------------------------------------------------%
    % Spacecraft state and parameters
    if t == 1
        r = r0;
        v = v0;
        qb = q0;
        wb = w0;
        mass_sc = Initial_SC_mass;
        time = t0;
    else
        r = r_tplus1;
        v = v_tplus1;
        qb = qb_tplus1;
        wb = wb_tplus1;
        time = tplus1;
    end

    Mass_Fuel  = Fuel_Mass - (1/(MR_LAM+1))*mass_spent_LAM - (1/(MR_ACT+1))*mass_spent_ACT;            % Current Mass of Fuel
    Mass_Oxi   = Oxi_Mass - (MR_LAM/(MR_LAM+1))*mass_spent_LAM - (MR_ACT/(MR_ACT+1))*mass_spent_ACT;   % Current Mass of Oxidizer


    % Store stuff at current timestep
    Output.States(t,:) = [(t-1)*tinc r' v' qb' wb'];
    Output.SC_Param(t,:) = [(t-1)*tinc mass_sc Mass_Fuel Mass_Oxi RCG_guess'];
    %---------------------------------------------------------------------%

    % LAM Disturbance
    F_LAM=[CartX;CartY;CartZ]*LAMThrust*LAM_Status;
    T_LAM = cross(-1*RCG_guess, F_LAM);

    % Reference attitude
    qref = qnorm(RefAttitude');

    % Controller (% Every 64mSec)
    %if(rem(time,0.064)<=1e-6 || t==1)
    if (rem((t-1), 4) == 0) || (t == 1)

        Qe=[qref(4) qref(3) -qref(2) -qref(1);
            -qref(3) qref(4) qref(1) -qref(2);
            qref(2) -qref(1) qref(4) -qref(3);
            qref(1) qref(2) qref(3) qref(4)]*QEstimate;

        Qe=qnorm(Qe);

        % Q error Saturation
        Qe(1) = max(min(Qe(1),qeLimit(1)),-qeLimit(1));
        Qe(2) = max(min(Qe(2),qeLimit(2)),-qeLimit(2));
        Qe(3) = max(min(Qe(3),qeLimit(3)),-qeLimit(3));

        we=(wbEstimate)-(wref);

        Tc=-Kp*[Qe(1);Qe(2);Qe(3)]-Kd*we;

        % Commanded Torque Saturation
        Tc(1) = max(min(Tc(1),TLimit(1)),-TLimit(1));
        Tc(2) = max(min(Tc(2),TLimit(2)),-TLimit(2));
        Tc(3) = max(min(Tc(3),TLimit(3)),-TLimit(3));

        % Normalization of the Commanded Torque
        C(1) = Tc(1)/TLimit(1);
        C(2) = Tc(2)/TLimit(2);
        C(3) = Tc(3)/TLimit(3);

    end

    %-----------------------DisturbanceTorque--------------------------------%
    if(rem(t,625)==0 || t==1)
        CurrentJD=JulianArray+((t-1)*tinc)/86400;
        et=cspice_str2et(sprintf('JD%0.9f',CurrentJD)); %Get ephemeris Time from SPICE function
        [psun]=cspice_spkpos('sun',et,'J2000','NONE','earth'); % Get Position of Sun from Earth in J2000 Frame
        SunVectorJ2000=psun/norm(psun);  % Unit Vector of Sun
        [a,e,i,AOP,Omega,TA]=Cartesian2Keplerian(r,v,'Eart');
        OrbitalElements = [e;a;i;Omega;AOP;TA];
        Text=dist_torque(CurrentJD,(t-1)*tinc,qb,IMsc,'Eart',OrbitalElements,RCG_guess,SunVectorJ2000,APanel,BPanel,CPanel);
    end


    % PWPFM (Every 16mSeconds)
    %if(rem(time,0.016)<=1e-6 || t==1)
    if (rem((t-1), 1) == 0) || (t == 1)
        for i=1:3
            y(i)=((C(i)-Um(i))*Km*tinc16)/Tm+py(i)*(1-tinc16/Tm);
            py(i)=y(i);
            if(y(i)>0)
                if(y(i)>Uon)
                    Um(i)=1;
                    LP(i)=1;
                elseif(y(i)<Uoff)
                    Um(i)=0;
                    LP(i)=0;
                else
                    Um(i)=LP(i);
                end
            end
            if(y(i)<0)
                if(y(i)<-Uon)
                    Um(i)=-1;
                    LNP(i)=-1;
                elseif(y(i)>-Uoff)
                    Um(i)=0;
                    LNP(i)=0;
                else
                    Um(i)=LNP(i);
                end
            end
        end
    end

    % TSL
    row=14-9*Um(1)-3*Um(2)-Um(3);
    MU=[T1(row) T4(row) T5(row) T8(row)];

    % TMatrix Calculuation
    % Torque and FMatrix = Direction cosines, due to each thruster, TotalThrusters = 4 in Normal mode
    TMatrix=zeros(3,TotalThrusters);
    FMatrix=zeros(3,TotalThrusters);
    for i=1:TotalThrusters
        Fx=cosd(CT(2,SelectedThrusters(i)))*cosd(CT(1,SelectedThrusters(i)));
        Fy=cosd(CT(2,SelectedThrusters(i)))*sind(CT(1,SelectedThrusters(i)));
        Fz=sind(CT(2,SelectedThrusters(i)));
        rx=RT(1,SelectedThrusters(i))-RCG_guess(1);                          % Thruster Position from Current CG Location
        ry=RT(2,SelectedThrusters(i))-RCG_guess(2);
        rz=RT(3,SelectedThrusters(i))-RCG_guess(3);
        TMatrix(:,i)=cross([rx;ry;rz],[Fx;Fy;Fz]);
        FMatrix(:,i)=[Fx;Fy;Fz];
    end

    % Net moment produced by ACT firing
    Sum_of_Thrusters = sum(MU);
    MUM = TMatrix*MU'*NominalThrust;

    % Net force produced by ACT firing
    ACTForce = FMatrix * MU' * NominalThrust;

    %TotalForce(:,1)=[F_LAM(1)+ACTForcex;F_LAM(2)+ACTForcey;F_LAM(3)+ACTForcez];
    TotalForce = F_LAM + ACTForce;
    BodyAcc = TotalForce/(mass_sc*1000);  % Inertial acceleration of center of mass due to thrust, in body frame, Km/s^2 (NOTE : May need to store this)
    R_ECI_BF = quat2dcm([qb(4),qb(1),qb(2),qb(3)]);
    Inertial_Acc = R_ECI_BF' * BodyAcc;   % Inertial acceleration of center of mass due to thrust, in inertial frame, km/s^2 (NOTE : May need to store this)

    % Position Velocity Propagation
    [y_tplus1] = RK4(tinc,time,[r;v],Inertial_Acc);
    r_tplus1 = y_tplus1(1:3);
    v_tplus1 = y_tplus1(4:6);

    % Attitude Propagation
    qb_tplus1 = q_dynamics(wb,qb,tinc);
    qb_tplus1 = qnorm(qb_tplus1);

    % Dynamics & Kinematics
    k11 = InvIMsc * (T_LAM + Text + MUM -(cross(wb,IMsc*wb)));
    H11 = wb + (tinc/2)*k11;
    k12 = InvIMsc * (T_LAM + Text + MUM -(cross(H11,IMsc*H11)));
    I11 = wb + 0.5*tinc*k12;
    k13 = InvIMsc * (T_LAM + Text + MUM -(cross(I11,IMsc*I11)));
    J11 = wb + tinc*k13;
    k14 = InvIMsc * (T_LAM + Text + MUM -(cross(J11,IMsc*J11)));

    wb_tplus1 = wb + (tinc/6)*(k11+2*k12+2*k13+k14);
    tplus1 = time + tinc

    %-----------------------GyroModel------------------------------%
    % Included on 17/12/15
    GyroRandomNumber=RandomNumbers(1:2,t);
    [WMeas,Bias]=GyroLandis(wb_tplus1,Bias,tinc,GyroRandomNumber);

%     w_lpf=gyro_lpf(WMeas,w_lpf,tinc,tau);
    wbEstimate=WMeas-BiasEstimate;

    % If SSU data is not available
    % Propagate using Gyro else use LKF
    if(rem(t,7)~=0)
        [QEstimate]=propagation(QEstimate,wbEstimate,tinc);
    else
        LKFRandomNumber=RandomNumbers(3:5,t);
        [QEstimate,BiasEstimate]=LKF(qb_tplus1,QEstimate,wbEstimate,BiasEstimate,tinc,Starnoise,KG,LKFRandomNumber);
    end

    % Update spacecraft parameters for next iteration
    mass_sc = mass_sc-(m_dot * tinc * Sum_of_Thrusters) - (m_dot_LAM*LAM_Status)*tinc; % This will be the mass @ the next time step
    mass_spent = mass_spent + (m_dot * tinc * Sum_of_Thrusters) + (m_dot_LAM)*tinc;
    mass_spent_LAM = mass_spent_LAM + (m_dot_LAM)*tinc;                 % Mass Spent by LAM from the beginning of the burn till next timestep
    mass_spent_ACT = mass_spent_ACT + (m_dot*tinc*Sum_of_Thrusters);    % Mass Spent by ACT from the beginning of the burn till next timestep

    Output.ACT_CMD(t,:) = [(t-1)*tinc MU]; % Store thruster commands
    Output.Control(t,:) = [(t-1)*tinc Qe']; % Store controller outputs

    % Post_Simulation_Calc
    qActual= qb_tplus1;
    qEstimate= QEstimate;
    Qe_Estimate =[ QEstimate(4) QEstimate(3) -QEstimate(2) -QEstimate(1);
        -QEstimate(3) QEstimate(4) QEstimate(1) -QEstimate(2);
        QEstimate(2) -QEstimate(1) QEstimate(4) -QEstimate(3);
        QEstimate(1) QEstimate(2) QEstimate(3) QEstimate(4)]*qActual;

    BiasEstimate_error = Bias-BiasEstimate;

    % Store estimate outputs
    Output.Estimate(t,:) = [(t-1)*tinc [Qe_Estimate(1)*114.8  Qe_Estimate(2)*114.8  Qe_Estimate(3)*114.8]*3600 BiasEstimate_error' BiasEstimate'];

end
Last_Step_Index = tfF;
cspice_kclear;
%% Plots
figure(1);
subplot(3,1,1); plot(Output.States(1:Last_Step_Index,1), Output.States(1:Last_Step_Index,12)*180/pi); ylabel('w_x (deg/s)');
subplot(3,1,2); plot(Output.States(1:Last_Step_Index,1), Output.States(1:Last_Step_Index,13)*180/pi); ylabel('w_y (deg/s)');
subplot(3,1,3); plot(Output.States(1:Last_Step_Index,1), Output.States(1:Last_Step_Index,14)*180/pi); ylabel('w_z (deg/s)');

figure(2);
subplot(3,1,1); plot(Output.States(1:Last_Step_Index,1), Output.States(1:Last_Step_Index,8 )); ylabel('q1 ');
subplot(3,1,2); plot(Output.States(1:Last_Step_Index,1), Output.States(1:Last_Step_Index,9 )); ylabel('q2 ');
subplot(3,1,3); plot(Output.States(1:Last_Step_Index,1), Output.States(1:Last_Step_Index,10)); ylabel('q3 ');

%{
figure(2);
[R1,R2,R3] = quat2angle([Output.Control(1:Last_Step_Index,5) Output.Control(1:Last_Step_Index,2) Output.Control(1:Last_Step_Index,2) Output.Control(1:Last_Step_Index,3)]);
subplot(3,1,1); plot(Output.States(1:Last_Step_Index,1), R1*180/pi); ylabel('Error Euler Angle Z');
subplot(3,1,2); plot(Output.States(1:Last_Step_Index,1), R2*180/pi); ylabel('Error Euler Angle Y');
subplot(3,1,3); plot(Output.States(1:Last_Step_Index,1), R3*180/pi); ylabel('Error Euler Angle X');

figure(3);
subplot(4,1,1); plot(Output.ACT_CMD(1:Last_Step_Index,1), Output.ACT_CMD(1:Last_Step_Index,2)); ylabel('T1');
subplot(4,1,2); plot(Output.ACT_CMD(1:Last_Step_Index,1), Output.ACT_CMD(1:Last_Step_Index,3)); ylabel('T4');
subplot(4,1,3); plot(Output.ACT_CMD(1:Last_Step_Index,1), Output.ACT_CMD(1:Last_Step_Index,4)); ylabel('T5');
subplot(4,1,4); plot(Output.ACT_CMD(1:Last_Step_Index,1), Output.ACT_CMD(1:Last_Step_Index,5)); ylabel('T8');

figure(4);
subplot(3,1,1);plot(Output.Estimate(1:Last_Step_Index,1)/1,Output.Estimate(1:Last_Step_Index,2),'m','linewidth',2);grid;zoom;legend('X-Axis')
title('Attitude Error','fontweight','b')
xlabel('Time,Seconds','fontweight','b')
ylabel('Arc Seconds','fontweight','b')
set(gca,'fontweight','b')
subplot(3,1,2);plot(Output.Estimate(1:Last_Step_Index,1)/1,Output.Estimate(1:Last_Step_Index,3),'m','linewidth',2);grid;zoom;legend('Y-Axis')
xlabel('Time,Seconds','fontweight','b')
ylabel('Arc Seconds','fontweight','b')
set(gca,'fontweight','b')
subplot(3,1,3);plot(Output.Estimate(1:Last_Step_Index,1)/1,Output.Estimate(1:Last_Step_Index,4),'m','linewidth',2);grid;zoom;legend('Z-Axis')
xlabel('Time,Seconds','fontweight','b');
ylabel('Arc Seconds','fontweight','b')
set(gca,'fontweight','b')

figure(5);
subplot(3,1,1); plot(Output.Estimate(1:Last_Step_Index,1)/3600,Output.Estimate(1:Last_Step_Index,5)*180/pi*3600,'m','linewidth',2);grid;zoom;legend('Yaw')
title('Error in Bias','fontweight','b')
xlabel('Time,Seconds','fontweight','b')
ylabel('Yaw Rate(deg/Hour)','fontweight','b')
set(gca,'fontweight','b')
subplot(3,1,2); plot(Output.Estimate(1:Last_Step_Index,1)/3600,Output.Estimate(1:Last_Step_Index,6)*180/pi*3600,'m','linewidth',2);grid;zoom;legend('Pitch')
xlabel('Time,Seconds','fontweight','b')
ylabel('Pitch Rate(deg/Hour)','fontweight','b')
set(gca,'fontweight','b')
subplot(3,1,3); plot(Output.Estimate(1:Last_Step_Index,1)/3600,Output.Estimate(1:Last_Step_Index,7)*180/pi*3600,'m','linewidth',2);grid;zoom;legend('Roll')
xlabel('Time,Seconds','fontweight','b')
ylabel('Roll Rate(deg/Hour)','fontweight','b')

figure(6);
subplot(3,1,1); plot(Output.Estimate(1:Last_Step_Index,1),Output.Estimate(1:Last_Step_Index,8)*180/pi,'m','linewidth',2);grid;zoom;legend('Yaw')
title('Bias','fontweight','b')
xlabel('Time,Seconds','fontweight','b')
ylabel('Yaw Rate(deg/Sec)','fontweight','b')
set(gca,'fontweight','b')
subplot(3,1,2); plot(Output.Estimate(1:Last_Step_Index,1),Output.Estimate(1:Last_Step_Index,9)*180/pi,'m','linewidth',2);grid;zoom;legend('Pitch')
xlabel('Time,Seconds','fontweight','b')
ylabel('Pitch Rate(deg/Sec)','fontweight','b')
set(gca,'fontweight','b')
subplot(3,1,3); plot(Output.Estimate(1:Last_Step_Index,1),Output.Estimate(1:Last_Step_Index,10)*180/pi,'m','linewidth',2);grid;zoom;legend('Roll')
xlabel('Time,Seconds','fontweight','b')
ylabel('Roll Rate(deg/Sec)','fontweight','b')

%}
