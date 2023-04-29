function qb=q_dynamics(wb,q,s)
F= 0.5*[0,  wb(3), -wb(2), wb(1);
        -wb(3), 0,  wb(1),  wb(2); 
         wb(2),-wb(1), 0,   wb(3);
        -wb(1), -wb(2), -wb(3), 0];
k1= F*q;
k2= F*(q+0.5*s*k1);
k3= F*(q+0.5*s*k2);
k4= F*(q+s*k3);
qbq = q+(s/6)*(k1+2*k2+2*k3+k4);
qb = qnorm(qbq);
end