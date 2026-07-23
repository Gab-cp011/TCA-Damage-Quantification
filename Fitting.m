function [q,err] = Fitting(H,f,f0,maxiter)

% compute_lcf_modalpar1 Compute modal parameters for LCF methods
% ------------------------------------------------------------------------------
% v = compute_lcf_modalpar1(H,f,f0,maxiter)
% ------------------------------------------------------------------------------
%
% Arguments : H  : measured FRF matrix of size nx1
%                  with n = nb of frequencies
%             f  : vector of frequencies (n-length)
%             f0 : vector of initial frequencies
%             maxiter : maximum iteration number
%
%             q  : vector of identified modal parameters
%                   
%
% -------------------------------------------- initialize
q = [];
err = 0;
if nargin ~= 4
  return
end
if isempty(f)
  return
end
n = length(H);
np = length(f0);  % nb of poles
w0 = 2*pi*f0(:);
w = 2*pi*f(:);
jw = sqrt(-1)*w;
% -------------------------------------------- initialize poles values
s0 = -w0/1000 + sqrt(-1)*w0;
jw0 = sqrt(-1) * mean(2*pi*f(:));
% -------------------------------------------- initialize modal parameter values
t0 = zeros(np,1);
q = [0 ;0 ; t0 ; s0];                  % ui,vi,ti,si
% -------------------------------------------- iterative loop
ir = 1;
iter = 1;
while (iter <= maxiter)
  % -------------------------------------------- build A matrix
  A = zeros(n,2+np);
  A(:,1) = -1;                                   % u
  A(:,2) = -(jw - jw0);                          % v
  for i=1:np
    A(:,2+i) = -1 ./ (jw - q(2+np+i));
  end
  
  if ir ~= 2
    for i=1:np
      A1 = -1 * (H -q(1) - q(2)*(jw - jw0));
      Atmp = zeros(n,np-1);
      for k=1:np
        if k ~= i 
          Atmp(:,k) = q(2+k) ./ (jw - q(2+np+k));
        end
      end
      A2 = sum(Atmp,2);
      A(:,2+np+i) = (A1 + A2) ./ (jw - q(2+np+i));
    end
  end
  % -------------------------------------------- build B matrix
  B1 = zeros(n,np);
  for i=1:np 
    B1(:,i) = q(2+i) ./ (jw - q(2+np+i));
  end
  B1 = sum(B1,2);
  B = -1 * (H -q(1) - q(2)*(jw - jw0)) + B1;
  % -------------------------------------------- norm matrices
  ncol = size(A,2);
  normA = sum(abs(A),1);
  for i=1:ncol
    A(:,i) = A(:,i)/(normA(i) + (normA(i)==0));
  end
  normB = sum(abs(B));
  B = B/normB;
  % -------------------------------------------- test rank
  r = rank(A);
  if r ~= ncol
    err = 1;
    return
  end
  % -------------------------------------------- solve system 
  dq = A\B;
  dq = dq./normA'*normB;
  % -------------------------------------------- update parameters
  q(1:ncol) = q(1:ncol) + dq;
  % -------------------------------------------- iteration increment
  iter = iter + 1;
  % -------------------------------------------- convergence criteria
  rns = sum(abs(dq./q(1:ncol)).^2);
  rns = rns/np;
  rnx = sum(abs(dq./q(1:ncol)));
  if rns <= 1e-5
    ir = 2;
  end
  if rnx <= 1e-5
    iter = iter+100*maxiter;
  end
  if iter >= maxiter
    err = 2;
  end
end
% -------------------------------------------- EOF

