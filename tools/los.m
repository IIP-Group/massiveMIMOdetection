% =========================================================================
% -- Part of "Data Detection in Massive MU-MIMO" Simulator
% -------------------------------------------------------------------------
% -- (c) 2020 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@ethz.ch and caoscar@ethz.ch
% =========================================================================

%% Line of sight channel generation.
%  Special thanks to Sven Jacobsson for this function
function H_pwm = los(par)

  U = par.MT;
  B = par.MR;

  c = 3e8; % speed of light [m/s]
  lambda =  c / 2e9; % carrier wavelength [m]
  delta = .5; % antenna spacing
    
  % place users randomly in a circular sector [-angSec/2,angSec/2],
  % using the Hash Slinging Slasher method
  % - generate angular separations of at least angSepMin
  UE_sep = zeros(U-1,1);
  par.angSec = 120;
  sectorAvail = par.angSec;
  par.angSepMin = 1;
  for uu = 1:U-1 % user loop
    UE_sep(uu) = unifrnd(par.angSepMin,sectorAvail/(U-uu));
    sectorAvail = sectorAvail - UE_sep(uu);
  end
  % - permute angular separations
  UE_sep = UE_sep(randperm(U-1));
  % - convert angular separations into angles
  UE_ang = [0; cumsum(UE_sep)];
  UE_ang = UE_ang - (max(UE_ang)-min(UE_ang))/2; % center users
  % -- wiggle users up to the angular space available
  UE_angLeft = par.angSec-(max(UE_ang)-min(UE_ang));
  UE_ang = UE_ang + unifrnd(-UE_angLeft/2,UE_angLeft/2);    
 
  % distance spread  
  par.d_spread = 0; % distance spread
  d_max = 150 + par.d_spread; % max user distance [m]
  d_min = 150 - par.d_spread; % min user distance [m]     
  if par.d_spread > 0
    d_avg = 2/3*(d_max^3-d_min^3)/(d_max^2-d_min^2); % avg user distance [m]
  else
    d_avg = d_max;
  end
  d_UE = sqrt((d_max^2-d_min^2)*rand(U,1) + d_min^2);% UE dist [m]
 
  broadside_BS_deg = 0; % broadside of BS antenna array  [deg]
  aod_UE = UE_ang + broadside_BS_deg;
 
  coord_BS = [0, 0]; % BS coord.
  coord_UE = ones(U,1)*coord_BS + ...
             (d_UE*ones(1,2)).*[cosd(aod_UE), sind(aod_UE)]; % UE coord.
     
  d_ant_BS = delta * lambda; % distance between BS antenna elements [m]
 
  % array rotation
  Omega_BS_deg = wrapTo360(90 - broadside_BS_deg); % BS array rotation [deg]
  Omega_BS_rad = pi/180 * Omega_BS_deg; % BS array rotation [rad]
 
  % antenna elem. coordinates
  x_BS = coord_BS(1) + d_ant_BS*((1:B) - (B+1)/2)*cos(pi-Omega_BS_rad);
  y_BS = coord_BS(2) + d_ant_BS*((1:B) - (B+1)/2)*sin(pi-Omega_BS_rad);
  x_UE = coord_UE(:,1);
  y_UE = coord_UE(:,2);
 
  % coordinates
  xx_BS = ones(U,1)*x_BS; yy_BS = ones(U,1)*y_BS;
  xx_UE = x_UE*ones(1,B); yy_UE = y_UE*ones(1,B);
     
  % reference distance
  d_ref = sqrt((xx_BS - xx_UE).^2 + (yy_BS - yy_UE).^2);
 
  % angles
  theta_BS = Omega_BS_rad - pi/2 + atan2((yy_UE-yy_BS),(xx_UE-xx_BS));
     
  % distances between UE and BS antenna elements
  dd_ant_BS = d_ant_BS*ones(U,1)*((1:B)-1); 
  tt_BS = theta_BS(:,1)*ones(1,B);
     
  % distance according to PWM model
  d_pwm = d_ref(:,1)*ones(1,B) - dd_ant_BS.*sin(tt_BS);
 
  % channel matrix
  H_pwm = d_avg./d_pwm .* exp(-1i*2*pi*d_pwm/lambda);
  H_pwm = H_pwm.';
    
end
 
function lon = wrapTo360(lon)
 
  positiveInput = (lon > 0);
  lon = mod(lon, 360);
  lon((lon == 0) & positiveInput) = 360;
 
end
