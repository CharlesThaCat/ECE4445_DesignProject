function Zmot = ZmotSim(f, fs, RE, QMS, QES)
%ZMOTSIM Summary of this function goes here
%   Detailed explanation goes here
RES = RE * QMS / QES;
Zmot = RES * (1/QMS)*(f*1j/fs) ./ ((f*1j/fs).*(f*1j/fs) + (1/QMS)*(f*1j/fs) + 1);

end

