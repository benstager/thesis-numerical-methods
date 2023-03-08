function h = stepsizeSelector(tspan, Nx, Nx_0, Nt_0, p, q)
%NUMSTEPSSELECTOR returns stepsize to balance spatial and
%temporal error. NT satisfies the equation
%
%   (Nx0/Nx)^p = (Nt0/Nt)^q
%
% where we assume that the temporal and spatial error are equal when Nx=Nx0
% and Nt=Nt0, and h = diff(tspan)/Nt
% Parameters:
%   tspan - 2x1 vector containing start and stop time
%   Nx - number of spatial points
%   Nx_0 - reference number of spatial points
%   Nt_0 - reference number of temporal points
%   p - spatial order
%   q - temporal order

h = diff(tspan) / numStepsSelector(Nx, Nx_0, Nt_0, p, q);

end