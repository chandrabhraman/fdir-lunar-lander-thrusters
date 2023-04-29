function q = qmult(q1, q2);

qmat = [ q2(4)  q2(3)   -q2(2)  q2(1);
         -q2(3) q2(4)   q2(1)   q2(2);
         q2(2)  -q2(1)  q2(4)   q2(3);
         -q2(1)  -q2(2) -q2(3)  q2(4)];

     
%      QMAT(I,I)=Q2(4)
%       QMAT(1,2)=Q2(3)
%       QMAT(1,3)=-Q2(2)
%       QMAT(2,3)=Q2(1)
%       QMAT(2,1)=-QMAT(1,2)
%       QMAT(3,1)=-QMAT(1,3)
%       QMAT(3,2)=-QMAT(2,3)
%       QMAT(1,4)=Q2(1)
%       QMAT(2,4)=Q2(2)
%       QMAT(3,4)=Q2(3)
%       QMAT(4,1)=-QMAT(1,4)
%       QMAT(4,2)=-QMAT(2,4)
%       QMAT(4,3)=-QMAT(3,4)


q = qmat*q1;