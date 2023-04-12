function [cs, sn] = givens_set(v1, v2)
%  if (v1 == 0)
%    cs = 0;
%    sn = 1;
%  else
    t = sqrt(v1^2 + v2^2);
%    cs = abs(v1) / t;
%    sn = cs * v2 / v1;
    cs = v1 / t;  % see http://www.netlib.org/eispack/comqr.f
    sn = v2 / t;
%  end
end

