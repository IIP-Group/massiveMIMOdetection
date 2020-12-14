% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% BEACHES denoiser
% Special thanks to Seyed Hadi Mirfashbafan for this function
% -- Seyed Hadi Mirfarshbafan, Alexandra Gallyas-Sanhueza, Ramina Ghods,
% -- and Christoph Studer, "Beamspace Channel Estimation for Massive MIMO
% -- mmWave Systems: Algorithms and VLSI Design," 
% -- IEEE Transactions on Circuits and Systems I: Regular Papers,
% -- 2020.
function Hest = BEACHES(par,Hn,N0)

  if par.quadriga==0
    error('You should not use BEACHES with an iid Rayleigh-fading channel.')
  end

  hdenoised = zeros(size(Hn));
  SURE = zeros(par.MR,1);  

  for uu=1:par.MT
      
    hnoisy = fft(Hn(:,uu))/sqrt(par.MR);
    N = length(hnoisy);
    
    % find optimal threshold for shrinkage function
    sorth = sort(abs(hnoisy),'ascend');
    cumsum = 0;
    cumsuminv = sum(1./abs(hnoisy));
    tau_opt = inf;
    suremin = inf;
    tau_low = 0;
    
    for bb = 1:N
      tau_high = sorth(bb);
      tau = max(tau_low,min(tau_high,N0/(2*(N-bb+1))*cumsuminv));
      tau_low = sorth(bb);
      SURE(bb) = cumsum + (N-bb+1)*tau^2 + N*N0 ...
        - 2*N0*(bb-1)-tau*N0*cumsuminv;
      cumsum = cumsum + sorth(bb).^2;
      cumsuminv = cumsuminv - 1/sorth(bb);
      if SURE(bb)<suremin
        suremin = SURE(bb);
        tau_opt = tau;
      end
    end
    
    hdenoised(:,uu) = (hnoisy./abs(hnoisy)).*max(abs(hnoisy)-tau_opt,0);
    
  end
  
  Hest = ifft(hdenoised)*sqrt(par.MR);

end
