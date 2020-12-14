% =========================================================================
% -- Data Detection in Massive MU-MIMO Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% -------------------------------------------------------------------------
% -- If you use this simulator or parts of it, then you must cite our
% -- papers:
%
% -- Oscar Castañeda, Tom Goldstein, and Christoph Studer,
% -- "Data Detection in Large Multi-Antenna Wireless Systems via
% -- Approximate Semidefinite Relaxation," 
% -- IEEE Transactions on Circuits and Systems I: Regular Papers,
% -- vol. 63, no. 12, pp. 2334-2346, Dec. 2016.
%
% -- Charles Jeon, Oscar Castañeda, and Christoph Studer
% -- "A 354 Mb/s 0.37 mm2 151 mW 32-User 256-QAM Near-MAP Soft-Input
% -- Soft-Output Massive MU-MIMO Data Detector in 28nm CMOS," 
% -- IEEE Solid-State Circuits Letters, vol. 2, no. 9, pp. 127-130, 
% -- Oct. 2019.
% =========================================================================

function detection_MIMO_sim(varargin)

  % -- set up default/custom parameters
  
  if isempty(varargin)
    
    disp('using default simulation settings and parameters...')
        
    % set default simulation parameters     
    par.runId = 0;         % simulation ID (used to reproduce results)
    par.MR = 32;           % receive antennas 
    par.MT = 16;           % transmit antennas (set not larger than MR!)         
    par.mod = 'QPSK';      % modulation type: 'BPSK','QPSK','16QAM','64QAM'             
    par.trials = 1e4;      % number of Monte-Carlo trials (transmissions)
    par.simName = ...      % simulation name (used for saving results)
      ['ERR_', num2str(par.MR), 'x', num2str(par.MT), '_', ...
        par.mod, '_', num2str(par.trials),'Trials'];
    par.SNRdB_list = ...   % list of SNR [dB] values to be simulated
      0:2:12;
    par.quadriga = 0;      % use quadriga channels
    par.los = 0;           % use line-of-sight (LoS) channel
    par.chest = 'PERF';    % channel estimator to use: Options:
                           % 'PERF', 'ML', 'BEACHES'
    par.detector = ...     % define detector(s) to be simulated. Options:
     {'SIMO','MMSE',...    % 'SIMO', 'ML', 'MRC', 'ZF', 'MMSE', 'SDR_RAND',
      'TASER_R','LAMA',... % 'SDR_R1', 'TASER', 'TASER_R', 'RBR', 'LAMA',
      'ADMIN',...          % 'ADMIN', 'BOX', 'OCD_MMSE', 'OCD_BOX', 
      'OCD_BOX','KBEST'};  % 'LR_LLL_DFE_rZF', 'KBEST'
                           % NOTE: 'ML', 'SDR_RAND', and 'SDR_R1' take a
                           %       long time if used for large systems 
                           % NOTE: 'SDR_RAND' and 'SDR_R1' requires CVX,
                           %       available here: 
                           %       http://cvxr.com/cvx/download/                         
                                   
    % SDR_RAND parameters -------------------------------------------------
    par.SDR_RAND.L = 50;   % Number of randomizations
                           
    % TASER parameters ----------------------------------------------------
    par.TASER.iters = 100;       % Number of TASER iterations
    par.TASER.alphaScale = 0.99; % Alpha scale for TASER's step size.
    %Step size used for different systems and iid Rayleigh:
    %-------------------------------------------------------------
    % MR / MT | Example system | Modulation | par.TASER.alphaScale
    % ratio   | (MRxMT)        | scheme     |
    %-------------------------------------------------------------
    % 1       | 32x32          | BPSK       | 0.99
    % 2       | 64x32          | BPSK       | 0.95
    % 4       | 128x32         | BPSK       | 0.8
    % 8       | 256x32         | BPSK       | 0.75
    % 1       | 32x32          | QPSK       | 0.99
    % 2       | 64x32          | QPSK       | 0.99
    % 4       | 128x32         | QPSK       | 0.99
    % 8       | 256x32         | QPSK       | 0.85
    %-------------------------------------------------------------
    %For LoS channels, you will need to tune par.TASER.alphaScale
    %For 32x16 LoS QPSK, 0.85 worked well
    
    % TASER_R parameters --------------------------------------------------
    par.TASER_R.iters = 100;       % Number of TASER iterations
    % Please check the comments in 'TASER parameters' for guidelines on
    % par.TASER_R.alphaScale
    par.TASER_R.alphaScale = 0.99; % Alpha scale for TASER's step size.    
    par.TASER_R.L = 50;            % Number of randomizations
    
    % RBR parameters ------------------------------------------------------
    par.RBR.iters = 20;       % Number of RBR iterations
    
    % LAMA parameters ---------------------------------------------------
    par.LAMA.iters = 30;         % Number of LAMA iterations
    par.LAMA.theta_tau_s = 0.5;  % Damping constant for moment variance param, (0,1)
    par.LAMA.theta_tau_z = 0.5;  % Damping constant for signal variance param, (0,1)
    
    % ADMIN parameters ----------------------------------------------------
    par.ADMIN.betaScale = 3;  % =1 returns biased MMSE on first iteration,
                              % but tuning may improve performance
    par.ADMIN.iters = 5;      % Number of ADMIN iterations                           
    par.ADMIN.gamma = 2;      % Step size (>0) for Lagrangian vector update
                              % A value <1 ensures convergence of ADMM, but
                              % larger values may improve performance
                              
    % BOX parameters ----------------------------------------------------
    par.BOX.iters = 10;       % Number of BOX iterations
    par.BOX.tau = 2^-7;       % Step size (>0) for gradient descent

    % OCD MMSE parameters -------------------------------------------------
    par.OCD_MMSE.iters = 10;   % Number of OCD_MMSE iterations      
    
    % OCD BOX parameters -------------------------------------------------
    par.OCD_BOX.iters = 10;    % Number of OCD_BOX iterations
    
    % K-BEST parameters ---------------------------------------------------
    par.KBEST.K = 5;          % Number of best nodes to consider at a time
    
  else
      
    disp('use custom simulation settings and parameters...')    
    par = varargin{1};     % only argument is par structure
    
  end

  % -- initialization
  
  addpath('./detectors')
  addpath('./tools')
  
  % load QuaDRiGa channels if they will be used
  if par.quadriga
    channelFile = '~/QuadrigaChannels/mmMAGIC_UMi_';
    if par.los
      channelFile = [channelFile 'LOS'];
    else
      channelFile = [channelFile 'NLOS'];
    end
    channelFile = [channelFile '_U' int2str(par.MT) '_B' int2str(par.MR)];
    channelFile = [channelFile '_W1-1_trials10000_Norm0p5to2_' ... 
                   int2str(par.runId)];
    loadedChannel = load(channelFile);
  end
  
  % use runId random seed (enables reproducibility)
  rng(par.runId,'twister'); 

  % set up Gray-mapped constellation alphabet (according to IEEE 802.11)
  switch (par.mod)
    case 'BPSK'
      par.symbols = [ -1 1 ];
    case 'QPSK' 
      par.symbols = [ -1-1i,-1+1i, ...
                      +1-1i,+1+1i ];
    case '16QAM'
      par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
                      -1-3i,-1-1i,-1+3i,-1+1i, ...
                      +3-3i,+3-1i,+3+3i,+3+1i, ...
                      +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM'
      par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
                      -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
                      -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
                      -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
                      +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
                      +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
                      +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
                      +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
                         
  end

  % extract average symbol energy
  par.Es = mean(abs(par.symbols).^2); 
  
  % precompute bit labels
  par.Q = log2(length(par.symbols)); % number of bits per symbol
  par.bits = de2bi(0:length(par.symbols)-1,par.Q,'left-msb');

  % track simulation time
  time_elapsed = 0;
  
  % -- start simulation 
  
  % initialize result arrays (detector x SNR)
  % vector error rate:
  res.VER = zeros(length(par.detector),length(par.SNRdB_list)); 
  % symbol error rate:
  res.SER = zeros(length(par.detector),length(par.SNRdB_list));
  % bit error rate:
  res.BER = zeros(length(par.detector),length(par.SNRdB_list));

  % generate random bit stream (antenna x bit x trial)
  bits = randi([0 1],par.MT,par.Q,par.trials);

  % trials loop
  tic
  for t=1:par.trials
  
    % generate transmit symbol
    idx = bi2de(bits(:,:,t),'left-msb')+1;
    s = par.symbols(idx).';
  
    % generate iid Gaussian channel matrix & noise vector
    n = sqrt(0.5)*(randn(par.MR,1)+1i*randn(par.MR,1));
    n_H = sqrt(0.5)*(randn(par.MR,par.MT)+1i*randn(par.MR,par.MT));
    if par.quadriga
      H = squeeze(loadedChannel.channels(t,:,:));
    else
      if par.los
        H = los(par); % we will use the planar wave model
      else
        H = sqrt(0.5)*(randn(par.MR,par.MT)+1i*randn(par.MR,par.MT));
      end
    end
    
    % transmit over noiseless channel (will be used later)
    x = H*s;
  
    % SNR loop
    for k=1:length(par.SNRdB_list)
      
      % compute noise variance 
      % (average SNR per receive antenna is: SNR=Es*||H||_F^2/MR/N0)
      N0 = par.Es*norm(H,'fro')^2*10^(-par.SNRdB_list(k)/10)/par.MR;
      
      % transmit data over noisy channel
      y = x+sqrt(N0)*n;
      
      % channel estimation
      switch (par.chest)
        case 'PERF'
          Hest = H;
        case 'ML'
          N0_chest = N0/par.Es/par.MT;
          Hest = H + sqrt(N0_chest)*n_H;
        case 'BEACHES'
          N0_chest = N0/par.Es/par.MT;
          Hn = H + sqrt(N0_chest)*n_H;
          Hest = BEACHES(par,Hn,N0_chest);
        otherwise
          error('par.detector type not defined.')
      end
    
      % algorithm loop      
      for d=1:length(par.detector)

        switch (par.detector{d})     % select algorithms
          case 'SIMO'                % SIMO lower bound detector
            [idxhat,bithat] = SIMO(par,Hest,y,s);
          case 'ML'                  % ML detection using sphere decoding
            [idxhat,bithat] = ML(par,Hest,y);
          case 'MRC'                 % unbiased MRC detection
            [idxhat,bithat] = MRC(par,Hest,y);
          case 'ZF'                  % unbiased ZF detection
            [idxhat,bithat] = ZF(par,Hest,y);
          case 'MMSE'                % unbiased MMSE detector
            [idxhat,bithat] = MMSE(par,Hest,y,N0);
          case 'SDR_RAND'            % detection via SDR with randomization
            srng = rng;
            [idxhat,bithat] = SDR_RAND(par,Hest,y);
            rng(srng);
          case 'SDR_R1'              % detection via SDR with rank-1 approx
            srng = rng;
            [idxhat,bithat] = SDR_R1(par,Hest,y);
            rng(srng);
          case 'TASER'               % TASER detector
            [idxhat,bithat] = TASER(par,Hest,y);
          case 'TASER_R'             % TASER detector
            srng = rng;
            [idxhat,bithat] = TASER_R(par,Hest,y);
            rng(srng);
          case 'RBR'                 % RBR detector
            [idxhat,bithat] = RBR(par,Hest,y);
          case 'LAMA'                % LAMA detector
            [idxhat,bithat] = LAMA(par,Hest,y,N0);
          case 'ADMIN'               % ADMIN detector
            [idxhat,bithat] = ADMIN(par,Hest,y,N0);
          case 'BOX'                 % BOX detector
            [idxhat,bithat] = BOX(par,Hest,y);
          case 'OCD_MMSE'            % OCD MMSE detector
            [idxhat,bithat] = OCD_MMSE(par,Hest,y,N0);
          case 'OCD_BOX'             % OCD BOX detector
            [idxhat,bithat] = OCD_BOX(par,Hest,y);
          case 'LR_LLL_DFE_rZF'      % LLL-LR-aided + DFE + rZF detector
            [idxhat,bithat] = LR_LLL_DFE_rZF(par,Hest,y);
          case 'KBEST'               % K-Best detector
            [idxhat,bithat] = KBEST(par,Hest,y);
          otherwise
            error('par.detector type not defined.')      
        end

        % -- compute error metrics
        err = (idx~=idxhat);
        res.VER(d,k) = res.VER(d,k) + any(err);
        res.SER(d,k) = res.SER(d,k) + sum(err)/par.MT;    
        res.BER(d,k) = res.BER(d,k) + ...
                         sum(sum(bits(:,:,t)~=bithat))/(par.MT*par.Q);                   
        
      end % algorithm loop
                 
    end % SNR loop    
    
    % keep track of simulation time    
    if toc>10
      time = toc;
      time_elapsed = time_elapsed + time;
      fprintf('estimated remaining simulation time: %3.0f min.\n', ...
                time_elapsed*(par.trials/t-1)/60);
      tic
    end      
  
  end % trials loop
  
  % normalize results
  res.VER = res.VER/par.trials;
  res.SER = res.SER/par.trials;
  res.BER = res.BER/par.trials;
  res.time_elapsed = time_elapsed;
  
  % -- save final results (par and res structures)

  save([ par.simName '_' num2str(par.runId) ],'par','res');
  
  % -- show results (generates fairly nice Matlab plot) 
      
  marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:'};
  figure(1)
  for d = 1:length(par.detector)
    if d==1
      semilogy(par.SNRdB_list,res.BER(d,:),marker_style{d},'LineWidth',2)
      hold on
    else
      semilogy(par.SNRdB_list,res.BER(d,:),marker_style{d},'LineWidth',2)
    end
  end
  hold off
  grid on
  xlabel('average SNR per receive antenna [dB]','FontSize',12)
  ylabel('uncoded bit error rate (BER)','FontSize',12)
  axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-4 1])
  legend(par.detector,'FontSize',12,'Interpreter','none')
  set(gca,'FontSize',12)
    
end
