function cdf = normcdf0 (x, m, s)
% NORMCDF  CDF of the normal distribution
%  CDF = normcdf(X, M, S) computes the cumulative distribution
%  function (CDF) at X of the normal distribution with mean M
%  and standard deviation S.
%
%  CDF = normcdf(X) is equivalent to CDF = normcdf(X, 0, 1)

% Adapted for Matlab (R) from GNU Octave 3.0.1
% Original file: statistics/distributions/normcdf.m
% Original author: TT <Teresa.Twaroch@ci.tuwien.ac.at>

% Copyright (C) 1995, 1996, 1997, 2005, 2006, 2007 Kurt Hornik
% Copyright (C) 2008-2009 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if (~ ((nargin == 1) || (nargin == 3)))
    error ('normcdf: you must give one or three arguments');
end

if (nargin == 1)
    m = 0;
    s = 1;
end

if (~isscalar (m) || ~isscalar (s))
    [retval, x, m, s] = common_size (x, m, s);
    if (retval > 0)
        error ('normcdf: x, m and s must be of common size or scalar');
    end
end

sz = size (x);
cdf = zeros (sz);

if (isscalar (m) && isscalar(s))
    if (find (isinf (m) | isnan (m) | ~(s >= 0) | ~(s < Inf)))
        cdf = NaN * ones (sz);
    else
        cdf =  stdnormal_cdf ((x - m) ./ s);
    end
else
    k = find (isinf (m) | isnan (m) | ~(s >= 0) | ~(s < Inf));
    if (any (k))
        cdf(k) = NaN;
    end

    k = find (~isinf (m) & ~isnan (m) & (s >= 0) & (s < Inf));
    if (any (k))
        cdf(k) = stdnormal_cdf ((x(k) - m(k)) ./ s(k));
    end
end

cdf((s == 0) & (x == m)) = 0.5;

end


function cdf = stdnormal_cdf (x)
% STDNORMAL_CDF  CDF of the standard normal distribution
%  CDF = stdnormal_cdf(X)
%  For each component of X, compute the CDF of the standard normal
%  distribution at X.

% Adapted for Matlab (R) from GNU Octave 3.0.1
% Original file: statistics/distributions/stdnormal_cdf.m
% Original author: KH <Kurt.Hornik@wu-wien.ac.at>

% Copyright (C) 1995, 1996, 1997, 1998, 2000, 2002, 2004, 2005, 2006,
%               2007 Kurt Hornik
% Copyright (C) 2008-2009 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if (nargin ~= 1)
    error('stdnormal_cdf: you should provide one argument');
end

sz = size (x);
if (numel(x) == 0)
    error ('stdnormal_cdf: x must not be empty');
end

cdf = (ones (sz) + erf (x / sqrt (2))) / 2;

end


function inv = norminv (x, m, s)
% NORMINV  Quantile function of the normal distribution
%  INV = norminv(X, M, S) computes the quantile (the inverse
%  of the CDF) at X of the normal distribution with mean M
%  and standard deviation S.
%
%  INV = norminv(X) is equivalent to INV = norminv(X, 0, 1)

% Adapted for Matlab (R) from GNU Octave 3.0.1
% Original file: statistics/distributions/norminv.m
% Original author: KH <Kurt.Hornik@wu-wien.ac.at>

% Copyright (C) 1995, 1996, 1997, 2005, 2006, 2007 Kurt Hornik
% Copyright (C) 2008-2009 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if (nargin ~= 1 && nargin ~= 3)
    error ('norminv: you must give one or three arguments');
end

if (nargin == 1)
    m = 0;
    s = 1;
end

if (~isscalar (m) || ~isscalar (s))
    [retval, x, m, s] = common_size (x, m, s);
    if (retval > 0)
        error ('norminv: x, m and s must be of common size or scalars');
    end
end

sz = size (x);
inv = zeros (sz);

if (isscalar (m) && isscalar (s))
    if (find (isinf (m) | isnan (m) | ~(s > 0) | ~(s < Inf)))
        inv = NaN * ones (sz);
    else
        inv =  m + s .* stdnormal_inv (x);
    end
else
    k = find (isinf (m) | isnan (m) | ~(s > 0) | ~(s < Inf));
    if (any (k))
        inv(k) = NaN;
    end

    k = find (~isinf (m) & ~isnan (m) & (s > 0) & (s < Inf));
    if (any (k))
        inv(k) = m(k) + s(k) .* stdnormal_inv (x(k));
    end
end

k = find ((s == 0) & (x > 0) & (x < 1));
if (any (k))
    inv(k) = m(k);
end

inv((s == 0) & (x == 0)) = -Inf;
inv((s == 0) & (x == 1)) = Inf;

end