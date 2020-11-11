function Z = LESim(f, LE, Le, n)
%LESIM Summary of this function goes here
%   Detailed explanation goes here

Z_Le = (2*pi*f*1j).^n * Le;
Z_LE = (2*pi*f*1j) * LE;

Z = (Z_Le .* Z_LE) ./ (Z_Le + Z_LE);
%Z = Z_Le;
end

