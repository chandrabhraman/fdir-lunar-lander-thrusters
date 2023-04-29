function [T] = store_in_window_dist(row, MU, disturbingAccelEstimate, disturbingAccel_i_j, T )
%{
 function: Stores (faultModeAcc - AccBfEstimate) in windows of corresponding fault modes

 Inputs :
 	row - Thruster pattern row used from TSL table
	MU - Thruster pattern fired
 	disturbingAccel_i_j - Accelerations corr. to fault modes - NominalAcc (linear, Angular)
 	disturbingAccelEstimate - Measured Accelerations - NominalAcc  , (linear, Angular)
 	T - Array of structures containing windows for each fault modes
%}
	for i = 1:length(T)

		% push the values for each fault modes
		currentWindowSize = size(T(i).window,2);
		if (MU(i) == 0)
			% If fault mode is not active
			T(i).window(:,currentWindowSize+1) = [ disturbingAccelEstimate; MU(i) ];
			T(i).active0(:,currentWindowSize+1) = [ zeros(6,1); MU(i) ];
		else
			% If fault mode is active
% 			disturbingAccel_i_j(:,row,i) - disturbingAccelEstimate
			T(i).window(:,currentWindowSize+1) = [ [disturbingAccel_i_j(:,row,i) - disturbingAccelEstimate]; MU(i) ];
			T(i).active0(:,currentWindowSize+1) = [ disturbingAccelEstimate; MU(i) ];
		end

		if (currentWindowSize >= T(i).windowSizeLimit)
			T(i).window(:,1) = []; % Pop the values if currentWindowSize exceeds windowSizeLimit
			T(i).active0(:,1) = [];
			ntimesActive = sum(T(i).window(7,:));
			ntimesinActive = T(i).windowSizeLimit - ntimesActive;

			% Calculate mean square error
			if ntimesActive == 0
				T(i).mse = 40;
			else
				% Active meanSquareError
				activeColumns = T(i).window(7,:) == 1;
				activeAcc = T(i).window(1:6, activeColumns );

				T(i).mseTotal(1) = sum(sum(activeAcc.^2))/ntimesActive;

				% Inactive meanSquareError
				inactiveColumns = T(i).window(7,:) == 0;
				inactiveAcc = T(i).window(1:6, inactiveColumns );

				T(i).mseTotal(2) = sum(sum(inactiveAcc.^2))/ntimesinActive;

				mseActive0 = sum(sum(T(i).active0(1:6,:).^2))/ntimesActive;
				% T(i).mse = T(i).mse/mseActive0;
				% if T(i).mse > 40
					% T(i).mse = 40;
				% end
			end
		end

		% if T(i).mse < T(i).threshold
			% T(i).count = T(i).count + 1;
			% if (T(2).mse < T(2).threshold)
			% 	pause
			% end
		% end
	end
end  % function
