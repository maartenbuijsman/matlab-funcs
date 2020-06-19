function [x] = gaussElim(A,b)
% File gaussElim.m
%   This subroutine will perform Gaussian elmination
%   on the matrix that you pass to it.
%   i.e., given A and b it can be used to find x,
%        Ax = b
%
%   To run this file you will need to specify several
%   things:
%   A - matrix for the left hand side.
%   b - vector for the right hand side
%
%   The routine will return the vector x.
%   ex: [x] = gaussElim(A,b)
%     this will perform Gaussian elminiation to find x.
%
%


 N = max(size(A));


 % Perform Gaussian Elimination

  for j=2:N,
       for i=j:N,
          m = A(i,j-1)/A(j-1,j-1);
          A(i,:) = A(i,:) - A(j-1,:)*m;
          b(i) = b(i) - m*b(j-1);
       end
  end

% Perform back substitution

x = zeros(N,1);
x(N) = b(N)/A(N,N);

for j=N-1:-1:1,
    x(j) = (b(j)-A(j,j+1:N)*x(j+1:N))/A(j,j);
end


