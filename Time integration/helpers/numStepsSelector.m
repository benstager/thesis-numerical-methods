function NT = numStepsSelector(Nx, Nx_0, Nt_0, p, q)
%NUMSTEPSSELECTOR returns number of timesteps to balance spatial and
%temporal error. NT satisfies the equation
%
%   (Nx0/Nx)^p = (Nt0/Nt)^q
%
% where we assume that the temporal and spatial error are equal when Nx=Nx0
% and Nt=Nt0.
% Parameters:
%   Nx - number of spatial points
%   Nx_0 - reference number of spatial points
%   Nt_0 - reference number of temporal points
%   p - spatial order
%   q - temporal order

NT = ceil(Nt_0 / (Nx_0 / Nx)^(p / q));

end