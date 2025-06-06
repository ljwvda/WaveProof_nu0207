function [wout,wfull] = conjugatesymmetric_ext_Wave(win,shape,symmetry)
% returns conjugate symmetric approximation of the input
% works for both intvals and floats

wfull = symmetrytofulltensor_ext_Wave(win,shape,symmetry);
wfull = (wfull+conj(wfull(end:-1:1,end:-1:1,end:-1:1,end:-1:1,:)))/2;
wout = fulltosymmetrytensor_ext(wfull,shape,symmetry);

end

