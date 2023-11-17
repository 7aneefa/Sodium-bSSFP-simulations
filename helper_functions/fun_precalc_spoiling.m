% function:     fun_precalc_spoiling
% purpose:      precalculate spoiling angles, used in 'fun_phase_spoil'
% inputs:   	number of isochromats
% outputs:      spoiling angles

% 01.03.2021 - f.kratzer@dkfz.de

%%
function [cos_spoil, sin_spoil] = fun_precalc_spoiling(num_of_spins)

    % equidistant spoiling angles for all isochromats
    phi = -pi:2*pi/num_of_spins:(pi-2*pi/num_of_spins); 
    
    cos_spoil = [ones(size(cos(1.*phi))); ones(size(cos(1.*phi))); ...
                cos(1.*phi); cos(1.*phi); ones(size(cos(1.*phi))); ...
                cos(1.*phi); cos(1.*phi); cos(2.*phi); cos(2.*phi); ...
                ones(size(cos(1.*phi))); cos(1.*phi); cos(1.*phi); ...
                cos(2.*phi); cos(2.*phi); cos(3.*phi); cos(3.*phi)];
    sin_spoil = [zeros(size(cos(1.*phi))); zeros(size(cos(1.*phi))); ...
                sin(1.*phi); sin(1.*phi); zeros(size(cos(1.*phi))); ...
                sin(1.*phi); sin(1.*phi); sin(2.*phi); sin(2.*phi); ...
                zeros(size(cos(1.*phi))); sin(1.*phi); sin(1.*phi); ...
                sin(2.*phi); sin(2.*phi); sin(3.*phi); sin(3.*phi)];