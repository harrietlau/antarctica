%% BAYESIAN FIT TO ANTARCTICA PENINSULA VALUES
%  HL Sept 2021
%  Fit Antarctic Peninsula to thermodynamic state.
%  Roughly follows "LAB_fitting_bayesian" project (or "EX" in comments)
%  in the VBR library

%  warning...: HL does not often code in MATLAB, apologies for bad MATLAB
%  code.

close all
clear all
clc

% Things to look out for... in "fit_seismic_observations_Lloyd", change
% period of band of observation to fit seismic band used.

% Run case in Fig. 7 of paper (i.e., 50 km, with whatever flow law):
% can change q_method to "xfit_mxw", "eburgers_psp", "xfit_premelt"
% in "outdir" make sure a directory "posterior" exists -- all posterior
% arrays will be output into this directory
outdir='AP/OUT_xfit_premelt_50km/';
q_method='xfit_premelt';
% stick with 50 km for this run
% these inputs were averaged external to this code
filenames.Vs = 'AP/IN/AP_Lloyd_Vs_50km.mat';
filenames.LAB = 'AP/IN/AP_LAB_50km.mat';
          
location.lat = -64.4759; %lat; 
location.lon = 299.5906 ; %lon; 
location.z_min = 200; % averaging min depth for asth.
location.z_max = 250; % averaging max depth for asth.
location.smooth_rad = 5;

fobs_LAB = 1.; % frequency to fit LAB model

posterior_A = fit_seismic_observations_Lloyd(filenames, location, q_method,outdir);
posterior_L = fit_plate_AN(filenames, location, q_method, posterior_A,fobs_LAB,outdir);


%% EXTRACT THERMODYNAMICAL PARAMETERS FOR PLATE MODEL

iplate = find(posterior_L.p_zPlate==max(posterior_L.p_zPlate));
iTp = find(posterior_L.p_Tp==max(posterior_L.p_Tp));
zplate = posterior_L.zPlate(iplate);
Tp = posterior_L.Tp(iTp);
% load all the models
load(posterior_L.Files.VBR_Box);

out_data = VBR(iTp,iplate);
clear VBR;

%% PRODUCE THE MECHANICAL PROPERTIES WITH FULL FREQUENCY RANGE

f = logspace(-15,0,75);
VBR.in = out_data.in;
VBR.in.SV.f = f;
[VBR] = VBR_spine(VBR);
fname = [outdir 'VBRout.mat'];
save(fname,'VBR');

%% SAVE VBR OUTPUT TO ASCII FILES

%  Thermodynamic parameters
ovar = VBR.in.SV.T_K;
save(strcat(outdir,'T_Tp',string(Tp),'.dat'),'ovar','-ascii');
ovar = VBR.in.SV.P_GPa;
save([outdir 'P.dat'],'ovar','-ascii');
ovar = VBR.in.SV.rho;
save([outdir 'rho.dat'],'ovar','-ascii');

%  Rheological parameters
ovar = VBR.in.SV.f;
save([outdir 'f.dat'],'ovar','-ascii');
ovar = VBR.out.elastic.anharmonic.Gu;
save([outdir 'Gu.dat'],'ovar','-ascii');
ovar = VBR.out.anelastic.xfit_mxw.J1;
save([outdir 'J1.dat'],'ovar','-ascii');
ovar = VBR.out.anelastic.xfit_mxw.J2;
save([outdir 'J2.dat'],'ovar','-ascii');
ovar = VBR.out.anelastic.xfit_mxw.Q;
save([outdir 'Q.dat'],'ovar','-ascii');
ovar = VBR.out.elastic.anharmonic.Vsu;
save([outdir 'vsu.dat'],'ovar','-ascii');
ovar = VBR.out.viscous.HK2003.diff.eta;
save([outdir 'eta0_diff_HK.dat'],'ovar','-ascii');


