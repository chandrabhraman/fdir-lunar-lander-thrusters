function [T, ntimesActive, ntimesinActive] = store_in_window(row, MU, AccBfEstimate6x1, faultModeAcc, NominalAccBf6x1, T )
%{
 function: Stores (faultModeAcc - AccBfEstimate) in windows of corresponding fault modes

 Inputs :
 	row - Thruster pattern row used from TSL table
	MU - Thruster pattern fired
	AccBfEstimate6x1 - [alpha; a]6x1 vector of measurement data
	faultModeAcc - [6x27x4] list containing 6x1 acceleration vector for each of the fault modes in deg/s^2 , m/s^2
	NominalAccBf6x1 - [alpha; a]6x1 vector of acceleration data Calculated from Dynamics
 	T - Array of structures containing windows for each fault modes
%}
	factor = 1000000; % Replace this factor by RCG_guess to keep them in same order

	ntimesActive = zeros(length(MU),1); % nThrustersx1
	ntimesinActive = zeros(length(MU),1);

	% Loop over all the fault modes
	for i = 1:length(T)

		%% push the values in windows for each fault modes to populate it
		currentWindowSize = size(T(i).window,2);
		if (MU(i) == 0)
			% If fault mode is inactive
			T(i).window(:,currentWindowSize+1) = [ [AccBfEstimate6x1 - NominalAccBf6x1]; MU(i) ];
			% T(i).active0(:,currentWindowSize+1) = [ zeros(6,1); MU(i) ];

		else
			% If fault mode is active
			T(i).window(:,currentWindowSize+1) = [ [AccBfEstimate6x1 - faultModeAcc(:,row,i)]; MU(i) ];
			len = size(T(i).debug,1);
			T(i).debug(len+1,2:7) = [AccBfEstimate6x1 - faultModeAcc(:,row,i)];
			% T(i).active0(:,currentWindowSize+1) = [ [AccBfEstimate6x1 - NominalAccBf6x1]; MU(i) ];

		end

		%% If windows is populated
		% Start calculation of MSE,
		% Start moving the window
		if (currentWindowSize >= T(i).windowSizeLimit)
			T(i).window(:,1) = []; % Pop the values if currentWindowSize exceeds windowSizeLimit

			% T(i).active0(:,1) = [];
			ntimesActive(i) = sum(T(i).window(7,:));
			ntimesinActive(i) = sum(T(i).window(7,:) == 0);

			% Calculate Active meanSquareError if Sufficient data is available for MSE calculation
			if ntimesActive(i) ~= 0

				activeColumns = T(i).window(7,:) == 1;
				activeAngularAcc = T(i).window(1:3, activeColumns );
				activeLinearAcc  = T(i).window(5:6, activeColumns );

				T(i).mse(:,1) = [ sum(sum( activeAngularAcc.^2 ));
				  				  sum(sum( activeLinearAcc.^2  ))]/ntimesActive(i);

				T(i).mse(2,1) = factor*T(i).mse(2,1);

				mseActive0angularAcc = sum(sum(T(i).active0(1:3,:).^2))/ntimesActive(i);
				mseActive0linearAcc  = sum(sum(T(i).active0(5:6,:).^2))/ntimesActive(i);

				% T(i).mse(1,1) = T(i).mse(1,1)/mseActive0angularAcc;
				% T(i).mse(2,1) = T(i).mse(2,1)/mseActive0linearAcc;

			end

			% Calculate InActive meanSquareError if Sufficient data is available for MSE calculation
			if ntimesinActive(i) ~= 0

				inactiveColumns = (T(i).window(7,:) == 0);
				inactiveAngularAcc = T(i).window(1:3, inactiveColumns );
				inactiveLinearAcc  = T(i).window(5:6, inactiveColumns );

				T(i).mse(:,2) = [ sum(sum( inactiveAngularAcc.^2 ));
				  				  sum(sum( inactiveLinearAcc.^2  ))]/ntimesinActive(i);

				T(i).mse(2,2) = factor*T(i).mse(2,2);


			end

			% Fault detection
			if (ntimesActive )
				multiplier = T(i).mse(:,1) < T(i).threshold(:,1);
				if (sum(multiplier) == 2 ) % fault detected
					T(i).count = T(i).count + multiplier*1;
					% Fault Exoneration
					multiplier = T(i).mse(:,2) > T(i).threshold(:,2);
					if (sum(multiplier) == 2)
						T(i).count = [0 0]';
					end
				end
			end % if ntimesActive

		end % if

	end % for
end  % function
