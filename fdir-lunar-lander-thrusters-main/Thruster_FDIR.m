% Normal Mode
% This code uses telemetry data as generated by 'SSU_Gyro_FDIR.m'
% So, first run 'SSU_Gyro_FDIR.m' as per desired Intensional thruster failure
% Introduce Intensional faults by putting 0's in MU (line no. 313 - 318 ) in the file SSU_Gyro_FDIR.m

clear all
close all
clc

%% Setup
addpath(genpath('Mice_data'));
addpath(genpath('OtherFunctions'));
addpath(genpath('FunctionsKinematicDynamics'));
addpath(genpath('DisturbanceTorqueModel'));
addpath(genpath('GyroModels'));
addpath(genpath('FDIR'));

% Load Mice Kernels
cspice_furnsh('Mice_data\MICE\kernels\lsk\naif0011.tls.pc')
cspice_furnsh('Mice_data\MICE\kernels\pck\pck00010.tpc')
cspice_furnsh('Mice_data\MICE\kernels\pck\earth_fixed.tf')
cspice_furnsh('Mice_data\MICE\kernels\pck\earth_720101_070426.bpc')
cspice_furnsh('Mice_data\MICE\kernels\spk\de421Planets.bsp')
cspice_furnsh('Mice_data\MICE\kernels\pck\moon_pa_de421_1900-2050.bpc')
cspice_furnsh('Mice_data\MICE\kernels\pck\moon_080317.tf')

load OrbitElements.mat
load Block1_4.mat
load ThrusterParameters.mat
load PWPFParameters.mat
load PositionCanting.mat
load FuelTankParameters.mat
load C2.mat
load telemetryReq.mat % as generated by SSU_Gyro_FDIR.m file

%% Simulation Parameters

tinc=0.016; % Time Step in Seconds
tfF=round(20/tinc);
tinc16 = 0.016;
t_latency_gyro = 0.007; % Seconds
t_latency_accel = 0.007; % Seconds
t_thrust_buildup = 0.025; % Seconds

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


% Initializing estimates
QEstimate = qnorm([0; 0; 0; 1]); %qnorm(RefAttitude');
wbEstimate = [0.1*pi/180; 0.1*pi/180; 0.1*pi/180]; % Angular velocity in bf
QEstimate_Gyro= qnorm([0; 0; 0; 1]); %qnorm(RefAttitude'); % Q estimate From Gyro
BiasEstimate=[0; 0; 0];
prevWbEstimate = wbEstimate;
angAccBfEstimate = [0; 0; 0]; % Measured Angular acceleration in bf

RandomNumbers=randn(5, tfF);

%% Spacecraft Parameters
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
               0.555904783422407]; %[m] Center of mass of spacecraft @ Launch

IMU_Accelerometer_SF_MA = [ 0.0001                  8.29e-05                  8.29e-05
                           8.29e-05                    0.0001                  8.29e-05
                           8.29e-05                  8.29e-05                    0.0001];
R_IMU = [-772.8
         -232.9032
          321.5598]*1e-3; % in [m] Position of IMU in body frame

R_body_to_IMU = [1 0 0
                 0 1 0
                 0 0 1]; % Orientation of IMU w.r.t body frame
R_IMU_to_body = R_body_to_IMU';


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
qeLimit =[0.1089 0.1414 0.1080]'; %'[1.34581881533101 1.72473867595819 1.31445993031359]';
%% Initial conditions

% Initial time
t0 = telemetryThrusterFDIR.BodyAccEstimate(1,1); %Time stamp of first acceleration Measurement

% Initial Position and Velocity
r0 = PosVel(1,1:3)';     % Position Vector in Inertial [Km]
v0 = PosVel(1,4:6)';     % Velocity Vector in Inertial [Km/s]

% Declare Initial values of BodyAcc & NominalAngularAccBf for plotting
NominalLinearAccBf = [0, 0, 0]';
NominalAngularAccBf = [0, 0, 0]';

% Initial Attitude
q0 = qnorm([0;0;0;1]); % qnorm(RefAttitude'); % Normalizing Quaternions

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
ntimesActive = zeros(1,TotalThrusters);
ntimesinActive = zeros(1,TotalThrusters);

py = [0 0 0]';
Um = [0 0 0]';
LP = [0 0 0]';
LNP = [0 0 0]';

%% Storage
Output.States   = zeros(tfF+1, 14); % timestamp, r', v', qb', wb'
Output.SC_Param = zeros(tfF+1, 7);
Output.ACT_CMD = zeros(tfF+1, 5);
Output.Control = zeros(tfF+1,5);


F_LAM=[CartX;CartY;CartZ]*LAMThrust*LAM_Status;

% Nominal acceleration 6x1  [alpha(deg/s^2), a (m/s^2)]
NominalAccelTSL = generate_nominal_accelerations('Block1_4.mat', RCG_guess, F_LAM, IMsc, NominalThrust, Text, Initial_SC_mass, w0);

% fault mode acceleration 6x1 [alpha (deg/s^2), a (m/s^2)]
faultModeAcc = generate_faultmode_accelerations('Block1_4.mat', RCG_guess, F_LAM, IMsc, NominalThrust, Text, Initial_SC_mass, w0);

% Thruster FDIR .. contains
ThrusterFDIR.NominalAccel = zeros(tfF+1, 7); % timestamp, NominalLinearAccBfel, NominalAngularAccBf
ThrusterFDIR.Measurements = zeros(tfF+1, 7); % timestamp, LinearAccBfEstimate, angAccBfEstimate
ThrusterFDIR.disturbingAccel_i_j = zeros(tfF+1,6,27,TotalThrusters); % fault_mode accelerations - Nominal accelerations
ThrusterFDIR.disturbingAccelEstimate = zeros(tfF+1, 7); % Measured acceleration - Nominal acceleration

% Populating TSL table
TSL = zeros(27,TotalThrusters); % Store TSL table of the block used in TSL
for i=1:TotalThrusters
	TSL(:,i) = eval(strcat('T',num2str(SelectedThrusters(i)))); % Store TSL table in T

	% Store mse for plotting
	ThrusterFDIR.mse(i).active = zeros(0,3); % t, alpha_mse, a_mse
	ThrusterFDIR.mse(i).inactive = zeros(0,3);

	% Create window for each of the fault modes
	T(i).window = zeros(7,0); % 6 for acceleration & 1 for active/inactive
	T(i).windowSizeLimit = 20;
	T(i).active0 = zeros(7,0);
	T(i).mse = [0, 0;  % Mean squre error [active angAcc,    Inactive angAcc;
				0, 0]; % 				  active linearAcc, Inactive linearAcc]
	T(i).threshold = [10 , 0.5;  % threshold [active angAcc,    Inactive angAcc;
				      5  , 10 ]; %            active linearAcc, Inactive linearAcc]
	T(i).count = [0, 0]'; % Thruster faults
	T(i).debug = zeros(0,7);
end



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

	TMatrix=zeros(3,TotalThrusters); % Torque produced due to each thruster
    FMatrix=zeros(3,TotalThrusters); % DC's of each thrust vectors
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

	% Load Um or MU
	row = telemetryThrusterFDIR.row(t,2);
	MU = TSL(row,:);

	% Net moment produced by ACT firing
    Sum_of_Thrusters = sum(MU);
    MUM = TMatrix*MU'*NominalThrust;

    % Net force produced by ACT firing
    ACTForce = FMatrix * MU' * NominalThrust;

    TotalForce = F_LAM + ACTForce;
    NominalLinearAccBf = TotalForce/(mass_sc);  % Inertial acceleration of center of mass due to thrust, in body frame, m/s^2
    R_ECI_BF = quat2dcm([qb(4),qb(1),qb(2),qb(3)]);
    Inertial_Acc = R_ECI_BF' * NominalLinearAccBf;   % Inertial acceleration of center of mass due to thrust, in inertial frame, m/s^2

    % Position Velocity Propagation
    [y_tplus1] = RK4(tinc,time,[r;v],Inertial_Acc);
    r_tplus1 = y_tplus1(1:3);
    v_tplus1 = y_tplus1(4:6);

    % Attitude Propagation
    qb_tplus1 = q_dynamics(wb,qb,tinc);
    qb_tplus1 = qnorm(qb_tplus1);

    % Dynamics & Kinematics
	NominalAngularAccBf = InvIMsc * (T_LAM + Text + MUM -(cross(wb,IMsc*wb)));
    k11 = NominalAngularAccBf;
    H11 = wb + (tinc/2)*k11;
    k12 = InvIMsc * (T_LAM + Text + MUM -(cross(H11,IMsc*H11)));
    I11 = wb + 0.5*tinc*k12;
    k13 = InvIMsc * (T_LAM + Text + MUM -(cross(I11,IMsc*I11)));
    J11 = wb + tinc*k13;
    k14 = InvIMsc * (T_LAM + Text + MUM -(cross(J11,IMsc*J11)));

    wb_tplus1 = wb + (tinc/6)*(k11+2*k12+2*k13+k14);
    tplus1 = time + tinc

	% Load Measured Accelerations & Angular velocities, Bias
	LinearAccBfEstimate = telemetryThrusterFDIR.BodyAccEstimate(t + (t_thrust_buildup + t_latency_accel)/tinc,2:4)'; % in m/s^2
	wbEstimate = telemetryThrusterFDIR.wbEstimate(t + (t_thrust_buildup + t_latency_gyro)/tinc,2:4)'; % Bias removed Measurement in rad/S^2
	angAccBfEstimate = (wbEstimate - prevWbEstimate)/tinc; % rad/s^2
	prevWbEstimate = wbEstimate;

	AccBfEstimate6x1 = [angAccBfEstimate*180/pi; LinearAccBfEstimate];   % deg/s^2 & m/s^2
	NominalAccBf6x1  = [NominalAngularAccBf*180/pi; NominalLinearAccBf]; % deg/s^2 & m/s^2

	% Deviation of measured acceleration from Nominal acceleration values
	disturbingAccelEstimate6x1 = AccBfEstimate6x1 - NominalAccBf6x1;
	mat6x27 = (ones(27,1)*NominalAccBf6x1')';

	% Deviation of faultModeAcc from Nominal acceleration values
	disturbingAccel_i_j = faultModeAcc - cat(3,mat6x27, mat6x27, mat6x27, mat6x27);

	% T = store_in_window_dist(row, MU, disturbingAccelEstimate6x1, disturbingAccel_i_j, T );
	for i = 1:TotalThrusters
		len(i) = size(T(i).debug,1);
    end

    % Populate window for each of the fault modes according to incoming
    % data
	[T, ntimesActive, ntimesinActive] = store_in_window(row, MU, AccBfEstimate6x1, faultModeAcc, NominalAccBf6x1, T);
	for i = 1:TotalThrusters
		if MU(i) == 1
			T(i).debug(len(i)+1,1) = t*tinc;
		end
	end

	ThrusterFDIR.NominalAccel(t,:) = [(t-1)*tinc NominalLinearAccBf' NominalAngularAccBf'];
	ThrusterFDIR.disturbingAccelEstimate(t,:) = [(t-1)*tinc disturbingAccelEstimate6x1'];
	ThrusterFDIR.disturbingAccel_i_j(t,:,:,:) = disturbingAccel_i_j;
	ThrusterFDIR.Measurements(t,:) = [(t-1)*tinc LinearAccBfEstimate' angAccBfEstimate'];
	Output.ACT_CMD(t,:) = [(t-1)*tinc MU];

	for i = 1:TotalThrusters
		if (t > T(i).windowSizeLimit )
			% Start storing after Sufficient data is stored in the windows

			if ntimesActive(i) > 0
				% If mse is available
				nrows = size(ThrusterFDIR.mse(i).active,1);
				ThrusterFDIR.mse(i).active(nrows+1,:)   = [t*tinc T(i).mse(1,1) T(i).mse(2,1)];
			end

			if ntimesinActive(i) > 0
				% If mse is available
				nrows = size(ThrusterFDIR.mse(i).inactive,1);
				ThrusterFDIR.mse(i).inactive(nrows+1,:) = [t*tinc T(i).mse(1,2) T(i).mse(2,2)];
			end
		end

		if (sum(T(i).count) >= 2 )
			fprintf('Thruster %d is faulty \n',i);
		end
	end


	% Update spacecraft parameters for next iteration
    mass_sc = mass_sc-(m_dot * tinc * Sum_of_Thrusters) - (m_dot_LAM*LAM_Status)*tinc; % This will be the mass @ the next time step
    mass_spent = mass_spent + (m_dot * tinc * Sum_of_Thrusters) + (m_dot_LAM)*tinc;
    mass_spent_LAM = mass_spent_LAM + (m_dot_LAM)*tinc;                 % Mass Spent by LAM from the beginning of the burn till next timestep
    mass_spent_ACT = mass_spent_ACT + (m_dot*tinc*Sum_of_Thrusters);    % Mass Spent by ACT from the beginning of the burn till next timestep

    % Post_Simulation_Calc
    qActual= qb_tplus1;
    qEstimate= QEstimate;
    Qe_Estimate =[ QEstimate(4)  QEstimate(3) -QEstimate(2) -QEstimate(1);
                  -QEstimate(3)  QEstimate(4)  QEstimate(1) -QEstimate(2);
                   QEstimate(2) -QEstimate(1)  QEstimate(4) -QEstimate(3);
                   QEstimate(1)  QEstimate(2)  QEstimate(3)  QEstimate(4)]*qActual;

	BiasEstimate_error = Bias-BiasEstimate;

	% Store estimate outputs
	Output.Estimate(t,:) = [(t-1)*tinc [Qe_Estimate(1)*114.8  Qe_Estimate(2)*114.8  Qe_Estimate(3)*114.8]*3600 BiasEstimate_error' BiasEstimate' LinearAccBfEstimate' angAccBfEstimate'];

end
Last_Step_Index = tfF;
cspice_kclear;

%% Plots

figure(1);
subplot(4,1,1); plot(Output.ACT_CMD(1:Last_Step_Index,1), Output.ACT_CMD(1:Last_Step_Index,2)); ylabel('T1');
subplot(4,1,2); plot(Output.ACT_CMD(1:Last_Step_Index,1), Output.ACT_CMD(1:Last_Step_Index,3)); ylabel('T4');
subplot(4,1,3); plot(Output.ACT_CMD(1:Last_Step_Index,1), Output.ACT_CMD(1:Last_Step_Index,4)); ylabel('T5');
subplot(4,1,4); plot(Output.ACT_CMD(1:Last_Step_Index,1), Output.ACT_CMD(1:Last_Step_Index,5)); ylabel('T8');

%%
figure(2);
subplot(4,1,1); plot(ThrusterFDIR.mse(1).active(:,1), ThrusterFDIR.mse(1).active(:,2),'--*r','linewidth',1,'markeredgecolor','k','markersize',2);grid on;legend('T1')
	title('Active, Mean square Error for alpha, window size = 20','fontweight','b')
	xlabel('Time,Seconds','fontweight','b');
	axis([0,Inf,0,Inf]);
subplot(4,1,2); plot(ThrusterFDIR.mse(2).active(:,1), ThrusterFDIR.mse(2).active(:,2),'--*r','linewidth',1,'markeredgecolor','k','markersize',2);grid on;legend('T4')
	xlabel('Time,Seconds','fontweight','b');
	axis([0,Inf,0,Inf]);
subplot(4,1,3); plot(ThrusterFDIR.mse(3).active(:,1), ThrusterFDIR.mse(3).active(:,2),'--*r','linewidth',1,'markeredgecolor','k','markersize',2);grid on;legend('T5')
	xlabel('Time,Seconds','fontweight','b');
	axis([0,Inf,0,Inf]);
subplot(4,1,4); plot(ThrusterFDIR.mse(4).active(:,1), ThrusterFDIR.mse(4).active(:,2),'--*r','linewidth',1,'markeredgecolor','k','markersize',2);grid on;legend('T8')
	xlabel('Time,Seconds','fontweight','b');
	axis([0,Inf,0,Inf]);

figure(3);
subplot(4,1,1); plot(ThrusterFDIR.mse(1).active(:,1), ThrusterFDIR.mse(1).active(:,3),'--*r','linewidth',1,'markeredgecolor','k','markersize',2);grid on;legend('T1')
	title('Active, Mean square Error for linear acc, window size = 20','fontweight','b')
	xlabel('Time,Seconds','fontweight','b');
	axis([0,Inf,0,Inf]);
subplot(4,1,2); plot(ThrusterFDIR.mse(2).active(:,1), ThrusterFDIR.mse(2).active(:,3),'--*r','linewidth',1,'markeredgecolor','k','markersize',2);grid on;legend('T4')
	xlabel('Time,Seconds','fontweight','b');
	axis([0,Inf,0,Inf]);
subplot(4,1,3); plot(ThrusterFDIR.mse(3).active(:,1), ThrusterFDIR.mse(3).active(:,3),'--*r','linewidth',1,'markeredgecolor','k','markersize',2);grid on;legend('T5')
	xlabel('Time,Seconds','fontweight','b');
	axis([0,Inf,0,Inf]);
subplot(4,1,4); plot(ThrusterFDIR.mse(4).active(:,1), ThrusterFDIR.mse(4).active(:,3),'--*r','linewidth',1,'markeredgecolor','k','markersize',2);grid on;legend('T8')
	xlabel('Time,Seconds','fontweight','b');
	axis([0,Inf,0,Inf]);


figure(4);
subplot(4,1,1); plot(ThrusterFDIR.mse(1).inactive(:,1), ThrusterFDIR.mse(1).inactive(:,2),'--og','linewidth',1,'markeredgecolor','k','markersize',2);grid on;legend('T1')
	title('InActive, Mean square Error for alpha, window size = 20','fontweight','b')
	xlabel('Time,Seconds','fontweight','b');
	axis([0,Inf,0,Inf]);
subplot(4,1,2); plot(ThrusterFDIR.mse(2).inactive(:,1), ThrusterFDIR.mse(2).inactive(:,2),'--og','linewidth',1,'markeredgecolor','k','markersize',2);grid on;legend('T4')
	xlabel('Time,Seconds','fontweight','b');
	axis([0,Inf,0,Inf]);
subplot(4,1,3); plot(ThrusterFDIR.mse(3).inactive(:,1), ThrusterFDIR.mse(3).inactive(:,2),'--og','linewidth',1,'markeredgecolor','k','markersize',2);grid on;legend('T5')
	xlabel('Time,Seconds','fontweight','b');
	axis([0,Inf,0,Inf]);
subplot(4,1,4); plot(ThrusterFDIR.mse(4).inactive(:,1), ThrusterFDIR.mse(4).inactive(:,2),'--og','linewidth',1,'markeredgecolor','k','markersize',2);grid on;legend('T8')
	xlabel('Time,Seconds','fontweight','b');
	axis([0,Inf,0,Inf]);

figure(5);
subplot(4,1,1); plot(ThrusterFDIR.mse(1).inactive(:,1), ThrusterFDIR.mse(1).inactive(:,3),'--og','linewidth',1,'markeredgecolor','k','markersize',2);grid on;legend('T1')
	title('InActive, Mean square Error for linear acc, window size = 20','fontweight','b')
	xlabel('Time,Seconds','fontweight','b');
	axis([0,Inf,0,Inf]);
subplot(4,1,2); plot(ThrusterFDIR.mse(2).inactive(:,1), ThrusterFDIR.mse(2).inactive(:,3),'--og','linewidth',1,'markeredgecolor','k','markersize',2);grid on;legend('T4')
	xlabel('Time,Seconds','fontweight','b');
	axis([0,Inf,0,Inf]);
subplot(4,1,3); plot(ThrusterFDIR.mse(3).inactive(:,1), ThrusterFDIR.mse(3).inactive(:,3),'--og','linewidth',1,'markeredgecolor','k','markersize',2);grid on;legend('T5')
	xlabel('Time,Seconds','fontweight','b');
	axis([0,Inf,0,Inf]);
subplot(4,1,4); plot(ThrusterFDIR.mse(4).inactive(:,1), ThrusterFDIR.mse(4).inactive(:,3),'--og','linewidth',1,'markeredgecolor','k','markersize',2);grid on;legend('T8')
	xlabel('Time,Seconds','fontweight','b');
	axis([0,Inf,0,Inf]);

%%
%{
figure(3);
subplot(3,1,1); plot(ThrusterFDIR.Measurements(1:Last_Step_Index,1), ThrusterFDIR.Measurements(1:Last_Step_Index,2));grid;zoom;legend('alpha_x')
	title('Measured Acceleration','fontweight','b');
	xlabel('Time,Seconds'); ylabel('m/s^2');
subplot(3,1,2); plot(ThrusterFDIR.Measurements(1:Last_Step_Index,1), ThrusterFDIR.Measurements(1:Last_Step_Index,3));grid;zoom;legend('alpha_x')
	xlabel('Time,Seconds'); ylabel('m/s^2');
subplot(3,1,3); plot(ThrusterFDIR.Measurements(1:Last_Step_Index,1), ThrusterFDIR.Measurements(1:Last_Step_Index,4));grid;zoom;legend('alpha_x')
	xlabel('Time,Seconds'); ylabel('m/s^2');

figure(4);
subplot(3,1,1); plot(ThrusterFDIR.Measurements(1:Last_Step_Index,1),ThrusterFDIR.Measurements(1:Last_Step_Index,5)*180/pi,'m','linewidth',2);grid;zoom;legend('alpha_x')
	title('Estimated Angular Acceleration','fontweight','b')
	xlabel('Time,Seconds','fontweight','b'); ylabel('deg/s^2')
subplot(3,1,2); plot(ThrusterFDIR.Measurements(1:Last_Step_Index,1),ThrusterFDIR.Measurements(1:Last_Step_Index,6)*180/pi,'m','linewidth',2);grid;zoom;legend('alpha_y')
	xlabel('Time,Seconds','fontweight','b'); ylabel('deg/s^2');
subplot(3,1,3); plot(ThrusterFDIR.Measurements(1:Last_Step_Index,1),ThrusterFDIR.Measurements(1:Last_Step_Index,7)*180/pi,'m','linewidth',2);grid;zoom;legend('alpha_z')
	xlabel('Time,Seconds','fontweight','b'); ylabel('deg/s^2');

figure(5);
subplot(3,1,1); plot(ThrusterFDIR.disturbingAccelEstimate(1:Last_Step_Index,1),ThrusterFDIR.disturbingAccelEstimate(1:Last_Step_Index,2),'m','linewidth',2);grid;zoom;legend('a_x')
	title('Estimated disturbing Acceleration','fontweight','b')
	xlabel('Time,Seconds','fontweight','b'); ylabel('m/s^2')
subplot(3,1,2); plot(ThrusterFDIR.disturbingAccelEstimate(1:Last_Step_Index,1),ThrusterFDIR.disturbingAccelEstimate(1:Last_Step_Index,3),'m','linewidth',2);grid;zoom;legend('a_y')
	xlabel('Time,Seconds','fontweight','b'); ylabel('m/s^2');
subplot(3,1,3); plot(ThrusterFDIR.disturbingAccelEstimate(1:Last_Step_Index,1),ThrusterFDIR.disturbingAccelEstimate(1:Last_Step_Index,4),'m','linewidth',2);grid;zoom;legend('a_z')
	xlabel('Time,Seconds','fontweight','b'); ylabel('m/s^2');

figure(6)
subplot(3,1,1); plot(ThrusterFDIR.disturbingAccelEstimate(1:Last_Step_Index,1),ThrusterFDIR.disturbingAccelEstimate(1:Last_Step_Index,5)*180/pi,'m','linewidth',2);grid;zoom;legend('alpha_x')
	title('Estimated Disturbing Angular Acceleration','fontweight','b')
	xlabel('Time,Seconds','fontweight','b'); ylabel('deg/s^2)')
subplot(3,1,2); plot(ThrusterFDIR.disturbingAccelEstimate(1:Last_Step_Index,1),ThrusterFDIR.disturbingAccelEstimate(1:Last_Step_Index,6)*180/pi,'m','linewidth',2);grid;zoom;legend('alpha_y')
	xlabel('Time,Seconds','fontweight','b'); ylabel('deg/s^2)');
subplot(3,1,3); plot(ThrusterFDIR.disturbingAccelEstimate(1:Last_Step_Index,1),ThrusterFDIR.disturbingAccelEstimate(1:Last_Step_Index,7)*180/pi,'m','linewidth',2);grid;zoom;legend('alpha_z')
	xlabel('Time,Seconds','fontweight','b'); ylabel('deg/s^2)');

figure(6);
subplot(3,1,1); plot(ThrusterFDIR.NominalAccel(1:Last_Step_Index,1), ThrusterFDIR.NominalAccel(1:Last_Step_Index,2), 'm','linewidth',2);grid;zoom;legend('a_x')
	title('Nominal acceleration','fontweight','b');
	xlabel('Time,Seconds'); ylabel('m/s^2');
subplot(3,1,2); plot(ThrusterFDIR.NominalAccel(1:Last_Step_Index,1), ThrusterFDIR.NominalAccel(1:Last_Step_Index,3), 'm','linewidth',2);grid;zoom;legend('a_y')
	xlabel('Time,Seconds'); ylabel('m/s^2');
subplot(3,1,3); plot(ThrusterFDIR.NominalAccel(1:Last_Step_Index,1), ThrusterFDIR.NominalAccel(1:Last_Step_Index,4), 'm','linewidth',2);grid;zoom;legend('a_z')
	xlabel('Time,Seconds'); ylabel('m/s^2');

figure(7);
subplot(3,1,1); plot(ThrusterFDIR.NominalAccel(1:Last_Step_Index,1), ThrusterFDIR.NominalAccel(1:Last_Step_Index,5), 'm','linewidth',2);grid;zoom;legend('alpha_x')
	title('Nominal Angular Acceleration','fontweight','b');
	xlabel('Time,Seconds'); ylabel('deg/s^2');
subplot(3,1,2); plot(ThrusterFDIR.NominalAccel(1:Last_Step_Index,1), ThrusterFDIR.NominalAccel(1:Last_Step_Index,6), 'm','linewidth',2);grid;zoom;legend('alpha_y')
	xlabel('Time,Seconds'); ylabel('deg/s^2');
subplot(3,1,3); plot(ThrusterFDIR.NominalAccel(1:Last_Step_Index,1), ThrusterFDIR.NominalAccel(1:Last_Step_Index,7), 'm','linewidth',2);grid;zoom;legend('alpha_z')
	xlabel('Time,Seconds'); ylabel('deg/s^2');
%}

%%
%{
% Plot disturbing accelerations
for i=1:4
	for j = 1:60
		figure(i)
		subplot(10,6,j); plot(Output.disturbingAccel_i_j(1:Last_Step_Index,1,:,:))

for i = 1:6
	subplot(6,1,i); plot(ThrusterFDIR.disturbingAccelEstimate(1:Last_Step_Index,1), ThrusterFDIR.disturbingAccel_i_j(1:Last_Step_Index,i,1,1));
	xlabel('Time,Seconds');
	if (i < 4)
		ylabel('m/s^2');
	else
		ylabel('deg/s^2');
	end
end


for i = 1:4
	figure(i+1)
	l = 1;
	for j = 1:10
		for k = 1:6
			subplot(10,6,l); plot(Output.disturbingAccelEstimate(1:Last_Step_Index,1), Output.disturbingAccel_i_j(1:Last_Step_Index,k,j,i));
			l = l+1;
		end
	end
end
%}
