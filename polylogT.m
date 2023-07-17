function y = polylogT(c,z,acc)
% POLYLOG - Computes the c-polylogarithm of z (Li_c)
%
% Usage:   y = polylog(c,z)
%          y = polylog(c,z,acc)
%
% Input:  |z|<1 : complex number defined on open unit disk
%          c    : base of polylogarithm
%          acc  : cutoff accuracy
%
% Output: y
%
% -------------------------------------------------------------------------
%  Copyright (C) 2009 Delft University of Technology
%    Faculty of Civil Engineering and Geosciences
%    Willem Ottevanger 
% 
%    Modified by Tobias Ambjornsson, Lund University for 
%    increased computational speed close to z=1
%    NOTE: for z>=0.999999 the asymptotic result is used
% -------------------------------------------------------------------------
if nargin == 2
   acc = eps;
end

y  = z;
y0 = y;
for j = 1:length(z);
   if z(j)<0.999999
      k = 1;
      err = 1;
      zk = z(j);
      while (abs(err)>acc);
         k = k+1;
         kn = k^c;
         zk = zk.*z(j);
         err = zk./kn;
         y(j) = y(j)+err;
      end
   else
       y(j) = zeta(c)+gamma(1-c)*(-log(z(j)))^(c-1);
   end;
end

