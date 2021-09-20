function predicted_vals = calc_LAB_AN(VbrBoxFile,settings,fobs)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% predicted_vals = calc_LAB(VbrBoxFile, settings) - using HL/BH zlab
%
% Calculates estimated LAB (from Q profile) for fixed zPlate and Tpot.
%
% Parameters:
% -----------
%      VbrBoxFile       filename for box with vbr results
%                          probability for each parameter combination
%
%       settings        structure settings for fitting, with the following
%                       required field names
%           set_Tp          return best zPlate for this potential temp.
%           q_method        Q method to use for fitting (must be in VBR 
%                           box), defaults to first available
%           per_bw_max      max period to use for fitting [s]
%                               default 30
%           per_bw_min      min period to use for fitting [s]
%                               default 10
%           dz_adi_km       depth range to average over to get asth Q [km]
%                               default 40
%           interp_pts      number of points to use for interpolation of 
%                           Z and Q_z
%                               default 1000
%           q_LAB_method    method for LAB search (see find_LAB_Q), 
%                               default 'Q_factor'
%           q_LAB_value     value to use in LAB search
%                               default 20 (i.e., Q_factor=20)
%
% Output:
% -------
%       predicted_vals  structure of predictions with the following fields
%           zLAB            seismic LAB depth [km]
%           meanVs          mean Vs in the depth range of interest [km/s]
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fprintf("\n\n\nFitting LAB and asthenosphere velocity\n")
  load(VbrBoxFile); % loads VBR structure
  settings=checkFittingSettings(settings,VBR); % process input args

  % predictions (from prior models):
  %   preds.zLAB_Q: LAB depth as measured by Q
  %   preds.meanVs average asthenosphere velocity
  
  
  predicted_vals =  predictObs(VBR, settings,fobs);
  
end

function settings = checkFittingSettings(settings,VBRbox)
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % define defaults, add to settings if field does not exist, calc freq range.
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % define defaults
  Defaults.per_bw_max=5.0;
  Defaults.per_bw_min=0.5;
  available_Q=fieldnames(VBRbox(1).out.anelastic);
  Defaults.q_method=available_Q{1};
  Defaults.dz_adi_km = 40;
  Defaults.interp_pts=1000;
  Defaults.q_LAB_method = 'Q_factor'; % method for LAB search (see find_LAB_Q)
  Defaults.q_LAB_value = 20; % value to use in LAB search

  % check input settings for fields, add defaults if they dont exist
  available_opts=fieldnames(Defaults);
  for i = 1:numel(available_opts)
    this_opt=available_opts{i};
    if ~isfield(this_opt,settings)
      settings.(this_opt)=Defaults.(this_opt);
    end
  end

  % calculate frequency range
  settings.freq_range=sort([1/settings.per_bw_max 1/settings.per_bw_min]);
end

function predictions =  predictObs(Box,settings,fobs)
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  predicts Z_LAB based on Q for a Box
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % allocate zLAB, meanVs for every Box
  Z_LAB_AN = zeros(size(Box));
  meanVs= zeros(size(Box));
%   meanQs= zeros(size(Box)); % ** Currently unused **

  % Get frequency axis (same for all boxes)
  freq=Box(1).in.SV.f;
  
  [~,ifreq]=min(abs(freq-fobs));  
  display(string(ifreq))
  
  % for interpolating Q, Vs
  nn_pts=settings.interp_pts;
  dz_adi_km=settings.dz_adi_km; % average Vs within dz_adi_km of LAB for adiabatic
  q_method=settings.q_method; % Q method for VBR
  q_LAB_method = settings.q_LAB_method; % method for LAB search
  q_LAB_value = settings.q_LAB_value; % value to use in LAB search
  
  for iBox = 1:numel(Box)
    Z_km = Box(iBox).Z_km;
    Qs_fz = Box(iBox).out.anelastic.(q_method).Q;
    Vs_fz = Box(iBox).out.anelastic.(q_method).V/1000; % m/s to km/s

    % 2D interplations of Q(z,f) and Vs(z,f) over depth and frequency
    % HL edit so we save freq and z interpolated points
    [ Vs_zf, freqi, zi ] = interp_FreqZ(Vs_fz,freq,nn_pts,...
                                         Z_km,nn_pts);
                                                         
    %% HL turned off:
    %Z_LAB_eta(iBox) = find_LAB_Q(Qs_z,Z_km_interp,'method',q_LAB_method, ...
    %                                    'value',q_LAB_value,'z_min_km',zplate);
    Z_LAB_AN(iBox) = find_LAB_AN(Box(iBox));
     
  end
  
  
  predictions.zLAB_Q=Z_LAB_AN;
  predictions.meanVs=meanVs;
  
end

function Z_LAB_AN = find_LAB_AN(Box)
    % For An, definition of LAB is same as ZPlate (thermally defined)
    Z_LAB_AN = Box.BoxParams.var2val;
end

    




