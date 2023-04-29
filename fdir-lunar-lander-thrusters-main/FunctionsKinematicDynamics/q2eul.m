function [ang]=q2eul(Q,IR1,IR2,IR3) % Ang in degree is the output
% DCM Verified -25/11/2014. 
% 3-2-1 - Shud be Verifed. Compare with Euler Angle,Quaternions and
% Transformation Matrices
dcm = zeros(3,3); 
      dcm(1,1) = Q(1)*Q(1)-Q(2)*Q(2)-Q(3)*Q(3)+ Q(4)*Q(4);
      dcm(1,2) = 2*(Q(1)*Q(2)+Q(3)*Q(4));
      dcm(1,3) = 2*(Q(1)*Q(3)-Q(2)*Q(4));
      dcm(2,1) = 2*(Q(1)*Q(2)-Q(3)*Q(4));
      dcm(2,2) = Q(2)*Q(2)-Q(1)*Q(1)-Q(3)*Q(3)+ Q(4)*Q(4);
      dcm(2,3) = 2*(Q(2)*Q(3)+Q(1)*Q(4));
      dcm(3,1) = 2*(Q(1)*Q(3)+Q(2)*Q(4));
      dcm(3,2) = 2*(Q(2)*Q(3)-Q(1)*Q(4));
      dcm(3,3) = Q(3)*Q(3)-Q(1)*Q(1)-Q(2)*Q(2)+Q(4)*Q(4) ;
      
      
 
if((IR1==1)&(IR2==2)&(IR3==3)) 
        temp = sqrt(dcm(3,2)*dcm(3,2)+dcm(3,3)*dcm(3,3));
        ang(1) =-atan2(dcm(3,2),dcm(3,3))  ;      
        ang(2) = atan2(dcm(3,1),temp);
        ang(3) =-atan2(dcm(2,1),dcm(1,1));
   elseif((IR1==1)&(IR2==3)&(IR3==2))  
        temp = sqrt(dcm(1,1)*dcm(1,1)+dcm(3,1)*dcm(3,1));
        ang(1) = atan2(dcm(2,3),dcm(2,2)) ;       
        ang(2) = atan2(dcm(3,1),dcm(1,1));
        ang(3) =-atan2(dcm(2,1),temp);
    elseif((IR1==2)&(IR2==3)&(IR3==1))  
        temp = sqrt(dcm(1,1)*dcm(1,1)+dcm(1,3)*dcm(1,3));
        ang(1) =-atan2(dcm(3,2),dcm(2,2));        
        ang(2) =-atan2(dcm(1,3),dcm(1,1));
        ang(3) = atan2(dcm(1,2),temp);
   elseif((IR1==2)&(IR2==1)&(IR3==3))  
        temp = sqrt(dcm(3,1)*dcm(3,1)+dcm(3,3)*dcm(3,3));
        ang(1) =-atan2(dcm(3,2),temp);        
        ang(2) = atan2(dcm(3,1),dcm(3,3));
        ang(3) = atan2(dcm(1,2),dcm(2,2));
   elseif((IR1==3)&(IR2==2)&(IR3==1))  
        temp = sqrt(dcm(2,3)*dcm(2,3)+dcm(3,3)*dcm(3,3));
        ang(1) = atan2(dcm(2,3),dcm(3,3));        
        ang(2) =-atan2(dcm(1,3),temp);
        ang(3) = atan2(dcm(1,2),dcm(1,1));
  elseif((IR1==3)&(IR2==1)&(IR3==2))   
        temp = sqrt(dcm(1,3)*dcm(1,3)+dcm(3,3)*dcm(3,3));
        ang(1) = atan2(dcm(2,3),temp) ;       
        ang(2) =-atan2(dcm(1,3),dcm(3,3));
        ang(3) =-atan2(dcm(2,1),dcm(2,2));
    end
    
       ang = ang' *180/pi;
      
      
 