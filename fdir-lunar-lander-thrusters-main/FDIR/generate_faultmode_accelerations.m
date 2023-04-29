function faultModeAcc = generate_faultmode_accelerations( block, RCG_guess, F_LAM, IMsc, NominalThrust, Text, mass_sc, wb)
%{
    Generates Nominal & Fault Mode acceleration values
	"This function uses equal Nominal thrust for each of the thrusters"
    Input :
		block - String containing name of thruster block used
        RCG_guess - Position vector of CG in body frame
		F_LAM - LAM disturbance forces
		IMsc - Ineria matrix
		NominalThrust - Nominal thrust of ACT
		Text - External disturbances forces
		mass_sc - mass of spacecraft
    Output :
        faultModeAcc - Acceleration for each of the possible fault modes

%}
	% Load the thruster block given and it's position, canting
    load(block)
	load PositionCanting.mat

    TSL = zeros(27,TotalThrusters);
    faultModeAcc = zeros(6, 27, TotalThrusters); % Acceleration [linear ;Angular] values for all faulty cases

	InvIMsc = inv(IMsc);
	T_LAM = cross(-1*RCG_guess, F_LAM);

    TMatrix=zeros(3,TotalThrusters); % Torque produced due to each thruster
    FMatrix=zeros(3,TotalThrusters); % DC's of each thrust vectors
    for i=1:TotalThrusters
        Fx=cosd(CT(2,SelectedThrusters(i)))*cosd(CT(1,SelectedThrusters(i)));
        Fy=cosd(CT(2,SelectedThrusters(i)))*sind(CT(1,SelectedThrusters(i)));
        Fz=sind(CT(2,SelectedThrusters(i)));
        rx=RT(1,SelectedThrusters(i))-RCG_guess(1);
        ry=RT(2,SelectedThrusters(i))-RCG_guess(2);
        rz=RT(3,SelectedThrusters(i))-RCG_guess(3);
        TMatrix(:,i)=cross([rx;ry;rz],[Fx;Fy;Fz]);
        FMatrix(:,i)=[Fx;Fy;Fz];
        TSL(:,i) = eval(strcat('T',num2str(SelectedThrusters(i)))); % Store TSL table in T
    end

    %--------------- loop over each of the fault modes --------------------%
    for ithruster = 1:TotalThrusters
        faultyRows = find(TSL(:,ithruster)); % Rows for which ith thruster is used

		%--------  loop to calculate faultModeAcc for each case of the fault mode ------ %
		for row = faultyRows'
			% row is the row no. for which ithruster is used in TSL
            MU = TSL(row,:);
			MU(ithruster) = 0;
			MUM = TMatrix*MU'*NominalThrust; % Net moment produced by ACT firing
            ACTForce = FMatrix * MU' * NominalThrust; % Net force produced by ACT firing
            TotalForce = F_LAM + ACTForce;

            LinearAccBf = TotalForce/(mass_sc);  % Inertial acceleration of center of mass due to thrust, in body frame,m/s^2
            AngularAccBf = InvIMsc * (T_LAM + Text + MUM -(cross(wb,IMsc*wb)));

			faultModeAcc(:, row, ithruster) = [AngularAccBf*180/pi; LinearAccBf];
        end
	end

end  % function
