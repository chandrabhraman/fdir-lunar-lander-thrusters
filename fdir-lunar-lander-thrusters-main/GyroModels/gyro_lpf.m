function w_f_n=gyro_lpf(w_meas,w_f,Ts,tau)
w_f_n=(Ts/(Ts+tau))*w_meas+(tau/(Ts+tau))*w_f;
end