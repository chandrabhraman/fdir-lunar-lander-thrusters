function [QEstimate]=propagation(QEstimatein,wbEstimate,s)
F= 0.5*[0,  wbEstimate(3), -wbEstimate(2), wbEstimate(1);
       -wbEstimate(3), 0,  wbEstimate(1),  wbEstimate(2);
        wbEstimate(2),-wbEstimate(1), 0,   wbEstimate(3);
       -wbEstimate(1), -wbEstimate(2), -wbEstimate(3), 0];

    q=QEstimatein;
    k1= F*q;
    k2= F*(q+0.5*s*k1);
    k3= F*(q+0.5*s*k2);
    k4= F*(q+s*k3);

    qbq = q+(s/6)*(k1+2*k2+2*k3+k4);
    QEstimate=qnorm(qbq);
