function [driverRadius,RE,VT,Vgen,frequency,magnitude,phase] = txtParser(path)
% function [A] = txtParser(path)
% txt file to local variable parser
% path: path of txt file

fileID = fopen(path,'r');

tline = fgetl(fileID);
% driver frame radius in meter
driverRadius = str2double(tline(1:2))/100; 
% voice coil resistance in ohm
RE = str2double(tline(19:22));
% test box volume in m^3
VT = str2double(tline(34))*0.0283168;

tline = fgetl(fileID);
% generator voltage in Vrms
Vgen = str2double(tline(8:14));

A = textscan(fileID, '%f %f %f','HeaderLines', 1);
% frequency in Hz
frequency = A{1,1};
% ZVC magnitude in ohm
magnitude = A{1,2};
% ZVC phase in degree
phase = A{1,3};

fclose(fileID);

end

