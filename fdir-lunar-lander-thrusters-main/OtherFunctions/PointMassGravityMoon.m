function Acc=PointMassGravityMoon(Pos,PositionEarthSun,PositionMoonEarth)

smu=132712440017.99;      % Gravitational Parameter of Sun
emu=398600.4415;      % Gravitational Parameter of Earth

        Sun2SC = Pos - PositionEarthSun;
        Earth2SC = Pos - PositionMoonEarth;

        vtmp1 = Pos - 2*PositionEarthSun;
        vtmp2 = Pos - 2*PositionMoonEarth;
        
        Store1 = Pos'*vtmp1;
        Store2 = PositionEarthSun'*PositionEarthSun;
        Store3 = Pos'*vtmp2;
        Store4 = PositionMoonEarth'*PositionMoonEarth;
        qsun = Store1/Store2;
        qearth = Store3/Store4;

        fsun = qsun * ((3 + 3* qsun + qsun * qsun)/ (1+ (1+ qsun)^1.5));
        fearth = qearth * ((3 + 3* qearth + qearth * qearth)/ (1+ (1+ qearth)^1.5));

        d3sun = norm(Sun2SC)^3;
        d3moon = norm(Earth2SC)^3;
        
        asun = -smu * (Pos + fsun * PositionEarthSun) / d3sun;
        aearth = -emu * (Pos + fearth * PositionMoonEarth) / d3moon;
   
        Acc=asun+aearth;
end

    