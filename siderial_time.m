function [GST, LST] = siderial_time(JD,lon)

%From Algorithm 2 of Vallado, p.62 of first edition
%some additional inspiration from https://smallsats.org/2013/04/14/local-sidereal-time/

J2000 = 2451545.0;
T_uti = (JD-J2000)/36525;

gst0 = 100.4606184+36000.77005361*T_uti+.00038793*T_uti^2-2.6E-8*T_uti^3; %Local Greenwich time at 00:00:00
gst0 = mod(gst0,360);

dt = datetime(JD,'ConvertFrom','juliandate');

GST = gst0 + 360.98564724*(hour(dt) +minute(dt)/60 + second(dt)/3600)/24; %Local Greenwich time at a different UTI other than 00:00:00
GST = mod(GST,360);
LST = GST+lon; %Local Siderial Time
end