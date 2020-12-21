% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% detection via the Row-By-Row (RBR) method
% -- Hoi-To Wai, Wing-Kin Ma, and Anthony Man-Cho So,
% -- "Cheap Semidefinite Relaxation MIMO Detection using Row-By-Row Block
% -- Coordinate Descent," 
% -- IEEE International Conference on Acoustics, Speech and Signal 
% -- Processing (ICASSP), May 2011, pp. 3256-3259.
function shat = RBR(par,H,y)

  % -- convert to real domain
  switch par.mod
    case 'QPSK'
      % -- convert to real domain
      yR = [ real(y) ; imag(y) ];
      HR = [ real(H) -imag(H) ; imag(H) real(H) ];   
      % -- preprocessing for SDR  
      C = [HR'*HR , -HR'*yR ; -yR'*HR yR'*yR ];
      N = 2*par.MT+1; 
    case 'BPSK'
      % -- convert to real domain
      yR = [ real(y) ; imag(y) ];
      HR = [ real(H) ; imag(H) ];  
      % -- preprocessing for SDR  
      C = [HR'*HR , -HR'*yR ; -yR'*HR yR'*yR ];
      N = par.MT+1; 
    otherwise
      error('modulation not supported')
  end

  % -- parameters
  sigma = 1e-2/N;  % Barrier parameter, with value as suggested in the
                   % paper. However, the authors' code uses
                   % sigma = 1e-2/(4*par.MT+1)
  
  % -- initialization
  X = eye(N);
  
  % -- RBR iterations
  for pp=1:par.RBR.iters
    for k=1:N      
      idxset = [1:k-1 k+1:N];      
      c = C(idxset,k);
      z = X(idxset,idxset)*c;      
      gamma = z'*c;      
      if gamma>0
        X(idxset,k) = -1/(2*gamma)*(sqrt(sigma^2+4*gamma)-sigma)*z;
      else
        X(idxset,k) = zeros(N-1,1);
      end
      % -- Ensure that X is symmetric (Very important!)
      X(k,k)=1;
      X(k,idxset) = X(idxset,k);      
    end
  end
  % -- Recover the data from the last column of X.
  shat = sign(X(1:N-1,N));
  
  switch par.mod
    case 'QPSK'
      shat = shat(1:par.MT,1)+1i*shat(par.MT+1:end,1);
    case 'BPSK'  
      shat = shat(1:par.MT,1);
    otherwise
      error('not supported')
  end  
  
end
