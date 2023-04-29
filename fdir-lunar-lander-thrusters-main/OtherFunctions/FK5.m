
% FK5 Reduction - Modified from dcmeci2ecef
function DCM=FK5(JD,deltaAT,deltaUT1,dNutation,polarMotion)
global NutationData
% UTC Contains - Year,Month,Day,Hour,Minute,Sec
% deltaUT1 - Difference Between UTC & UTC1.
% detlaAT - Difference Between UTC & Atomic Time. 
% dNutation - dPsi, dEpsilon.
% polarMotion - X,Y

jdTT = JD + (deltaAT+32.184)/86400;
% Number of Julian centuries since J2000 for terrestrial time.
tTT = (jdTT - 2451545)/36525;
tTT2 = tTT*tTT;
tTT3 = tTT2*tTT;
tTT4 = tTT3*tTT;
% Additional time calculations:
        jdUT1 = JD + deltaUT1/86400;
        tUT1 = (jdUT1 - 2451545)/36525;
        tUT12 = tUT1*tUT1;
        tUT13 = tUT12*tUT1;
       
        ThetaGMST = (1.00965822615e6 + 4.746600277219299e10*tUT1 + 1.396560*tUT12 + 9.3e-5*tUT13)/3600;   % Degress
        thGMST = mod(ThetaGMST,360);

        
% Nutation
    % Mean obliquity of the ecliptic
        epsilonBar = (23.439291 - 0.0130042*tTT - 1.64e-7*tTT2 + 5.04e-7*tTT3)*pi/180;
        
        MAM=(134.96340251 + (1717915923.2178*tTT+31.8792*tTT2 + 0.051635*tTT3 - 0.00024470*tTT4)/3600)*pi/180;  % Mean Anamoly of Moon
        MAS =(357.52910918 + (129596581.0481*tTT - 0.5532*tTT2 - 0.000136*tTT3 - 0.00001149*tTT4)/3600)*pi/180; % Mean Anamoly of Sun.
        MArM = (93.27209062 + (1739527262.8478*tTT - 12.7512*tTT2 + 0.001037*tTT3 + 0.00000417*tTT4)/3600)*pi/180; % Mean Argument of Moon
        DLSM =(297.85019547 + (1602961601.2090*tTT - 6.3706*tTT2 + 0.006593*tTT3 - 0.00003169*tTT4)/3600)*pi/180; % Difference b/w Mean Longitude of Sun And Moon
        AMO =(125.04455501 + (-6962890.2665*tTT + 7.4722*tTT2 + 0.007702*tTT3 - 0.00005939*tTT4)/3600)*pi/180;  % Ascending Node of Moon Orbit
        
        NLong =0;   NObli=0;
        for i =1:106
            Temp1 = NutationData(i,1)*MAM + NutationData(i,2)*MAS + NutationData(i,3)*MArM +NutationData(i,4)*DLSM + NutationData(i,5)*AMO;
            Temp2 = (NutationData(i,6) + NutationData(i,7)*tTT)*sin(Temp1);
            Temp3 = (NutationData(i,8) + NutationData(i,9)*tTT)*cos(Temp1);
            NLong = NLong + Temp2; %Nutation in Longtitude 
            NObli = NObli + Temp3; %Nutation in Obliquity
        end
        nutationAngles=[NLong NObli];
      % Nutation angles obtained using JPL data
      
% This Function tries to Load De405 and gets Nutation Angle from it. Slow. 
%         nutationAngles = earthNutation(2400000.5+jdTT);       
        dpsi = nutationAngles(1); %Nutation in Longitude
        depsilon = nutationAngles(2); %Nutation in obliquity
        %The last two terms for equation of the equinoxes are only included
        %if the date is later than January 1, 1997 (MJD=50449)
        omegaMoon = (125.04455501 - (5*360 + 134.1361851)*tTT + 0.0020756*tTT2 + 2.139e-6*tTT3)*pi/180;
        omegaMoon(jdUT1<50449) = 0;
        % Equation of the equinoxes
        equinoxEq = dpsi.*cos(epsilonBar) + (0.00264/3600)*(pi/180)*sin(omegaMoon) + (0.000063/3600)*(pi/180)*sin(2*omegaMoon); 
        % Greenwhich apparent sidereal time
        thGAST = thGMST*pi/180 + equinoxEq;
        % Transformation matrix for earth rotation
        R = angle2dcm(thGAST,0,0);
        % Adjustments to nutation angles (provided from real measurements)
        dpsi = dpsi + dNutation(1);
        depsilon = depsilon + dNutation(2);
        % True obliquity of ecliptic
        epsilon = epsilonBar+depsilon;
        % True equator to mean equinox date transformation matrix
        N = angle2dcm(epsilonBar,-dpsi,-epsilon,'XZX');
% Precession
        % Zeta, theta and z represent the combined effects of general
        % precession.
        zeta = ((2306.2181*tTT + 0.30188*tTT2 + 0.017998*tTT3)/3600)*pi/180; 
        theta = ((2004.3109*tTT - 0.42665*tTT2 - 0.041833*tTT3)/3600)*pi/180;
        z = ((2306.2181*tTT + 1.09468*tTT2 + 0.018203*tTT3)/3600)*pi/180; 
        % Mean equinox to celestial reference frame
        P = angle2dcm(-zeta,theta,-z,'ZYZ');
        
  % Polar motion
        W = eye(3);
        W(1,3) =  polarMotion(1);
        W(3,1) = -polarMotion(1);
        W(2,3) = -polarMotion(2);
        W(3,2) =  polarMotion(2);
        
DCM=W*R*N*P;

end

% UnUsed
        
% Sidereal time - Was In  the Matlab Function dcm2ecef - Not Sure about
% jdUT1 and few calculations inside about the source, replaced it with the
% GMAT Equations. Works the same after 24 Hours same Results are obtained. 

%         % Greenwich mean sidereal time at midnight
%         thGMST0h = 100.4606184 + 36000.77005361*tUT1 + 0.00038793*tUT12 - 2.6e-8*tUT13;
%         % Ratio of universal to sidereal time
%         omegaPrec = 1.002737909350795 + 5.9006e-11*tUT1 - 5.9e-15*tUT12;
%         % Elapsed universal time since midnight to the time of the
%         % observation
%         UT1 = hour*60*60 + min*60 + ssUT;
%         % Greenwich mean sidereal time at time of the observation
%         thGMST = mod(thGMST0h + (360/(24*60*60))*omegaPrec.*UT1,360);


%- Replaced UTC times with Julian Date
% year=UTC(1);
% month = UTC(2);
% day = UTC(3);
% hour = UTC(4);
% min = UTC(5);
% sec = UTC(6);

% % Seconds for UT1
% ssUT = sec + deltaUT1;
% % Seconds for UTC
% ssTT = sec + deltaAT + 32.184;
% % Julian date for terrestrial time
% jdTT = mjuliandate(year,month,day,hour,min,ssTT);
% jdUT1 = mjuliandate(year,month,day,hour,min,ssUT);