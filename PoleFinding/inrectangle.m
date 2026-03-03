function [val,inds] = inrectangle(val, xmin,xmax,ymin,ymax)
    inds = real(val) > xmin & real(val) < xmax & imag(val) > ymin & imag(val) < ymax;
    val = val(inds);
end