function [val] = lapdygreen(x,y,ny)
% in two dimension, the derivitive of the Green's function in the ny
% direction.
squrr = (x(1)-y(1))^2+(x(2)-y(2))^2;
val = (x-y) * ny' /(2*pi*squrr);
end

