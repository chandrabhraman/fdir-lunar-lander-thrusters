% JD to Gregodrian Date

function [h,min,s]=JD2UTC(Julian)

T1900=(Julian-2415019.5)/365.25;
Year=1900+floor(T1900);
LeapYrs=floor((Year-1900-1)*(0.25));
Days=(Julian-2415019.5)-((Year-1900)(365)+LeapYrs);

if Days<1
    Year=Year-1;
    LeapYrs=floor((Year-1900-1)*(0.25));
    Days=(Julian-2415019.5)-((Year-1900)*(365)+LeapYrs);
end
if(mod(Year,4)==0)
    LMonth(2)=29;
end
DayofYr=floor(Days);

    