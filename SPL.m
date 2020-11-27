function spl = SPL(prms)
%SPL Calculate the SPL from the given prms value

spl = 20 * log10(prms / (2*10^(-5)));

end