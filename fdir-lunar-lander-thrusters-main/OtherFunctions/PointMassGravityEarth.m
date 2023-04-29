function Acc=PointMassGravityEarth(Pos,PositionEarthSun,PositionEarthMoon)

smu=132712440017.99;      % Gravitational Parameter of Sun
mmu=4902.8005821478;      % Gravitational Parameter of Moon

% PositionEarthSun=planetEphemeris(JD,'Earth','Sun');
% PositionEarthMoon=planetEphemeris(JD,'Earth','Moon');
% [~,~,PositionEarthSun]=SunPos(JD);
% PositionEarthSun=PositionEarthSun';
% [~,~,PositionEarthMoon]=MoonPos(JD);
% PositionEarthMoon=PositionEarthMoon';

% - Unused ---
%           SC2Sun=PositionEarthSun - Pos;
%           SC2Moon=PositionEarthMoon - Pos;
%           
%          asun = smu*(SC2Sun/norm(SC2Sun)^3 - PositionEarthSun/norm(PositionEarthSun)^3);
%          amoon = mmu*(SC2Moon/norm(SC2Moon)^3); %- PositionEarthMoon/norm(PositionEarthMoon)^3);
%---
        Sun2SC = Pos - PositionEarthSun;
        Moon2SC = Pos - PositionEarthMoon;

        vtmp1 = Pos - 2*PositionEarthSun;
        vtmp2 = Pos - 2*PositionEarthMoon;
        
        Store1 = Pos'*vtmp1;
        Store2 = PositionEarthSun'*PositionEarthSun;
        Store3 = Pos'*vtmp2;
        Store4 = PositionEarthMoon'*PositionEarthMoon;
        qsun = Store1/Store2;
        qmoon = Store3/Store4;

        fsun = qsun * ((3 + 3* qsun + qsun * qsun)/ (1+ (1+ qsun)^1.5));
        fmoon = qmoon * ((3 + 3* qmoon + qmoon * qmoon)/ (1+ (1+ qmoon)^1.5));

        d3sun = norm(Sun2SC)^3;
        d3moon = norm(Moon2SC)^3;
        
        asun = -smu * (Pos + fsun * PositionEarthSun) / d3sun;
        amoon = -mmu * (Pos + fmoon * PositionEarthMoon) / d3moon;
   
        Acc=asun+amoon;
end

function [R]=SunPos1(jdate)

TUT1 = (jdate-2451545)/36525;
LamdaM = (280.4606184+36000.77005361*TUT1)*pi/180;
M = (357.5277233 + 35999.05034*TUT1)*pi/180;
LamdaE= LamdaM + (1.914666471*pi/180)*sin(M) + 0.019994643*sin(2*M);
Epsilon = 23.439291*pi/180 - 0.0130042*TUT1;
R = [cos(LamdaE);cos(Epsilon)*sin(LamdaE);sin(Epsilon)*sin(LamdaE)];
end

function [rasc, decl, rsun] = SunPos (jdate)

% Solar ephemeris

% Input    
%  jdate = julian day

% output
%  rasc = right ascension of the sun (radians)
%         (0 <= rasc <= 2 pi)
%  decl = declination of the sun (radians)
%         (-pi/2 <= decl <= pi/2)
%  rsun = eci position vector of the sun (kilometers)

atr = pi / 648000;

rsun = zeros(1, 3);

% time arguments

djd = jdate - 2451545;     

t = (djd / 36525) + 1;

% fundamental arguments (radians)

gs = r2r(0.993126 + 0.0027377785 * djd);
lm = r2r(0.606434 + 0.03660110129 * djd);
ls = r2r(0.779072 + 0.00273790931 * djd);
g2 = r2r(0.140023 + 0.00445036173 * djd);
g4 = r2r(0.053856 + 0.00145561327 * djd);
g5 = r2r(0.056531 + 0.00023080893 * djd);
rm = r2r(0.347343 - 0.00014709391 * djd);

% geocentric, ecliptic longitude of the sun (radians)

plon = 6910 * sin(gs) + 72 * sin(2 * gs) - 17 * t * sin(gs);
plon = plon - 7 * cos(gs - g5) + 6 * sin(lm - ls) ... 
       + 5 * sin(4 * gs - 8 * g4 + 3 * g5);
plon = plon - 5 * cos(2 * (gs - g2)) - 4 * (sin(gs - g2) ... 
       - cos(4 * gs - 8 * g4 + 3 * g5));
plon = plon + 3 * (sin(2 * (gs - g2)) - sin(g5) - sin(2 * (gs - g5)));
plon = ls + atr * (plon - 17 * sin(rm));

% geocentric distance of the sun (kilometers)

rsm = 149597870.691 * (1.00014 - 0.01675 * cos(gs) - 0.00014 * cos(2 * gs));

% obliquity of the ecliptic (radians)

obliq = atr * (84428 - 47 * t + 9 * cos(rm));

% geocentric, equatorial right ascension and declination (radians)

a = sin(plon) * cos(obliq);
b = cos(plon);
   
rasc = atan3(a, b);
decl = asin(sin(obliq) * sin(plon));

% geocentric position vector of the sun (kilometers)

rsun(1) = rsm * cos(rasc) * cos(decl);
rsun(2) = rsm * sin(rasc) * cos(decl);
rsun(3) = rsm * sin(decl);

end

function [rasc, decl, rmoon] = MoonPos (jdate)

% lunar ephemeris

% input
%  jdate = julian date

% output
%  rasc  = right ascension of the moon (radians)
%          (0 <= rasc <= 2 pi)
%  decl  = declination of the moon (radians)
%          (-pi/2 <= decl <= pi/2)
%  rmoon = eci position vector of the moon (kilometers)

atr = pi / 648000;

rmoon = zeros(1, 3);

% time arguments

djd = jdate - 2451545;

t = (djd / 36525) + 1;

% fundamental arguments (radians)

gm = r2r(0.374897 + 0.03629164709 * djd);
gm2 = 2 * gm;
gm3 = 3 * gm;
fm = r2r(0.259091 + 0.0367481952 * djd);
fm2 = 2 * fm;
em = r2r(0.827362 + 0.03386319198 * djd);
em2 = 2 * em;
em4 = 4 * em;
gs = r2r(0.993126 + 0.0027377785 * djd);
lv = r2r(0.505498 + 0.00445046867 * djd);
lm = r2r(0.606434 + 0.03660110129 * djd);
ls = r2r(0.779072 + 0.00273790931 * djd);
rm = r2r(0.347343 - 0.00014709391 * djd);

% geocentric, ecliptic longitude of the moon (radians)

l = 22640 * sin(gm) - 4586 * sin(gm - em2) + 2370 * sin(em2);
l = l + 769 * sin(gm2) - 668 * sin(gs) - 412 * sin(fm2);
l = l - 212 * sin(gm2 - em2) - 206 * sin(gm - em2 + gs);
l = l + 192 * sin(gm + em2) + 165 * sin(em2 - gs);
l = l + 148 * sin(gm - gs) - 125 * sin(em) - 110 * sin(gm + gs);
l = l - 55 * sin(fm2 - em2) - 45 * sin(gm + fm2) + 40 * sin(gm - fm2);
l = l - 38 * sin(gm - em4) + 36 * sin(gm3) - 31 * sin(gm2 - em4);
l = l + 28 * sin(gm - em2 - gs) - 24 * sin(em2 + gs) + 19 * sin(gm - em);
l = l + 18 * sin(em + gs) + 15 * sin(gm + em2 - gs) + 14 * sin(gm2 + em2);
l = l + 14 * sin(em4) - 13 * sin(gm3 - em2) - 17 * sin(rm);
l = l - 11 * sin(gm + 16 * ls - 18 * lv) + 10 * sin(gm2 - gs) ... 
    + 9 * sin(gm - fm2 - em2);
l = l + 9 * (cos(gm + 16 * ls - 18 * lv) - sin(gm2 - em2 + gs)) ... 
    - 8 * sin(gm + em);
l = l + 8 * (sin(2 * (em - gs)) - sin(gm2 + gs)) - 7 * (sin(2 * gs) ... 
    + sin(gm - 2 * (em - gs)) - sin(rm));
l = l - 6 * (sin(gm - fm2 + em2) + sin(fm2 + em2)) ...
    - 4 * (sin(gm - em4 + gs) - t * cos(gm + 16 * ls - 18 * lv));
l = l - 4 * (sin(gm2 + fm2) - t * sin(gm + 16 * ls - 18 * lv));
l = l + 3 * (sin(gm - 3 * em) - sin(gm + em2 + gs) ... 
    - sin(gm2 - em4 + gs) + sin(gm - 2 * gs) + sin(gm - em2 - 2 * gs));
l = l - 2 * (sin(gm2 - em2 - gs) + sin(fm2 - em2 + gs) - sin(gm + em4));
l = l + 2 * (sin(4 * gm) + sin(em4 - gs) + sin(gm2 - em));

plon = lm + atr * l;

% geocentric, ecliptic latitude of the moon (radians)

b = 18461 * sin(fm) + 1010 * sin(gm + fm) + 1000 * sin(gm - fm);
b = b - 624 * sin(fm - em2) - 199 * sin(gm - fm - em2) ... 
    - 167 * sin(gm + fm - em2);
b = b + 117 * sin(fm + em2) + 62 * sin(gm2 + fm) + 33 * sin(gm - fm + em2);
b = b + 32 * sin(gm2 - fm) - 30 * sin(fm - em2 + gs) ... 
    - 16 * sin(gm2 - em2 + fm);
b = b + 15 * sin(gm + fm + em2) + 12 * sin(fm - em2 - gs) ... 
    - 9 * sin(gm - fm - em2 + gs);
b = b - 8 * (sin(fm + rm) - sin(fm + em2 - gs)) ... 
    - 7 * sin(gm + fm - em2 + gs);
b = b + 7 * (sin(gm + fm - gs) - sin(gm + fm - em4));
b = b - 6 * (sin(fm + gs) + sin(3 * fm) - sin(gm - fm - gs));
b = b - 5 * (sin(fm + em) + sin(gm + fm + gs) + sin(gm - fm + gs) ... 
    - sin(fm - gs) - sin(fm - em));
b = b + 4 * (sin(gm3 + fm) - sin(fm - em4)) - 3 * (sin(gm - fm - em4) ... 
    - sin(gm - 3 * fm));
b = b - 2 * (sin(gm2 - fm - em4) + sin(3 * fm - em2) - sin(gm2 - fm + em2) ... 
    - sin(gm - fm + em2 - gs));

plat = atr * (b + 2 * (sin(gm2 - fm - em2) + sin(gm3 - fm)));

% obliquity of the ecliptic (radians)

obliq = atr * (84428 - 47 * t + 9 * cos(rm));

% geocentric distance (kilometers)

r = 60.36298 - 3.27746 * cos(gm) - .57994 * cos(gm - em2);
r = r - .46357 * cos(em2) - .08904 * cos(gm2) + .03865 * cos(gm2 - em2);
r = r - .03237 * cos(em2 - gs) - .02688 * cos(gm + em2) ... 
    - .02358 * cos(gm - em2 + gs);
r = r - .0203 * cos(gm - gs) + .01719 * cos(em) + .01671 * cos(gm + gs);
r = r + .01247 * cos(gm - fm2) + .00704 * cos(gs) + .00529 * cos(em2 + gs);
r = r - .00524 * cos(gm - em4) + .00398 * cos(gm - em2 - gs) ... 
    - .00366 * cos(gm3);
r = r - .00295 * cos(gm2 - em4) - .00263 * cos(em + gs) ... 
    + .00249 * cos(gm3 - em2);
r = r - .00221 * cos(gm + em2 - gs) + .00185 * cos(fm2 - em2) ... 
    - .00161 * cos(2 * (em - gs));
r = r + 0.00147 * cos(gm + fm2 - em2) - 0.00142 * cos(em4) ... 
    + 0.00139 * cos(gm2 - em2 + gs);

rmm = 6378.14 * (r - 0.00118 * cos(gm - em4 + gs) - 0.00116 * cos(gm2 + em2) ... 
      - 0.0011 * cos(gm2 - gs));

% geocentric, equatorial right ascension and declination (radians)

a = sin(plon) * cos(obliq) - tan(plat) * sin(obliq);
b = cos(plon);

rasc = atan3(a, b);

decl = asin(sin(plat) * cos(obliq) + cos(plat) * sin(obliq) * sin(plon));

% geocentric position vector of the moon (kilometers)

rmoon(1) = rmm * cos(rasc) * cos(decl);
rmoon(2) = rmm * sin(rasc) * cos(decl);
rmoon(3) = rmm * sin(decl);

end
    
    