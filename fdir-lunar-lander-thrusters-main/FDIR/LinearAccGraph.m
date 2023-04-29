for i = 1:4
	activeRows = (TSL(:,i) == 1);
	a = 1:28;
	figure()
	plot(a(activeRows), faultModeAcc(4,activeRows,i));
	title('acceleration_x'); ylabel('m/s^2');
	figure()
	plot(a(activeRows), faultModeAcc(5,activeRows,i));
	title('acceleration_y'); ylabel('m/s^2');
	figure()
	plot(a(activeRows), faultModeAcc(6,activeRows,i));
	title('acceleration_z'); ylabel('m/s^2');
	pause
end
