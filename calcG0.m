function [x] = calcG0(w)
%calcG0 calculates G0(0,w)
%   given w, calculates G0(0,w)



t = 1.0;
eta = 1e-3; %width of delta peaks
w = w+1i*eta;

x = (t./(sqrt(w-2*t).*sqrt(w+2*t)));


end

