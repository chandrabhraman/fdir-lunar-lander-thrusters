function [year, month, day, hour, minute, second] = jd2jdate(jday)
% Modified to Give output as UTCG
%JD2JDATE Julian day number to UTCG Date.
%
%   [YEAR, MONTH, DAY, HOUR, MINUTE, SECOND] = JD2JDATE(JDAY) returns the
%   Gregorian calendar date (year, month, day, hour, minute, and second)
%   
%
%   Start of the JD (Julian day) count is from 0 at 12 noon 1 JAN -4712
%   (4713 BC), Julian proleptic calendar.  Note that this day count conforms
%   with the astronomical convention starting the day at noon, in contrast
%   with the civil practice where the day starts with midnight.
%
%   Astronomers have used the Julian period to assign a unique number to
%   every day since 1 January 4713 BC.  This is the so-called Julian Day
%   (JD). JD 0 designates the 24 hours from noon UTC on 1 January 4713 BC
%   (Julian calendar) to noon UTC on 2 January 4713 BC.

%   Sources:  - http://tycho.usno.navy.mil/mjd.html
%             - http://www.tondering.dk/claus/cal/julperiod.php



   nargsin = nargin;
   error(nargchk(1, 1, nargsin));

   ijday = floor(jday);                 % integer part
   fjday = jday - ijday;                % fraction part

   second = 86400 * fjday;
   hour   = floor(second / 3600);         % get number of hours
   second = second - 3600 * hour;         % remove the hours
   minute = floor(second / 60);           % get number of minutes
   second = floor(second - 60 * minute);         % remove the minutes
   hour   = hour + 12;                  % Julian days start at noon

   % The following algorithm is from the Calendar FAQ.
    a = ijday+32044;
    b = floor((4*a+3)/(146097));
    c = a-floor(146097*b/4);


   d = floor((4 * c + 3) / 1461);
   e = c - floor((1461 * d) / 4);
   m = floor((5 * e + 2) / 153);

   day   = e - floor((153 * m + 2) / 5) + 1;
   month = m + 3 - 12 * floor(m / 10);
   year  = b * 100 + d - 4800 + floor(m / 10);