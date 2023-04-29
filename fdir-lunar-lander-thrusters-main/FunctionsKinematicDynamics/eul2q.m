function [q] = eul2q(ang,IR1,IR2,IR3)   % ang in deg is input

q1 =[0;0;0;1];
q2 =[0;0;0;1];
q3 =[0;0;0;1];      
q1(1) = sin(ang(1)*pi/360.0);
q1(4) = cos(ang(1)*pi/360.0);

q2(2) = sin(ang(2)*pi/360.0);
q2(4) = cos(ang(2)*pi/360.0);

q3(3) = sin(ang(3)*pi/360.0);
q3(4) = cos(ang(3)*pi/360.0);
if((IR1==1)&&(IR2==2)&&(IR3==3))
[qtemp] = qmult(q1,q2);
[q] = qmult(qtemp,q3);
elseif((IR1==1)&&(IR2==3)&&(IR3==2))
[qtemp] = qmult(q1,q3);
[q] = qmult(qtemp,q2);
elseif((IR1==2)&&(IR2==3)&&(IR3==1))
[qtemp] = qmult(q2,q3);
[q] = qmult(qtemp,q1);
elseif((IR1==2)&&(IR2==1)&&(IR3==3))
[qtemp] = qmult(q2,q1);
[q] = qmult(qtemp,q3);
elseif((IR1==3)&&(IR2==2)&&(IR3==1))
[qtemp] = qmult(q3,q2);
[q] = qmult(qtemp,q1);
elseif((IR1==3)&&(IR2==1)&&(IR3==2))
[qtemp] = qmult(q3,q1);
[q] = qmult(qtemp,q2);
end
[q] = qnorm(q);
