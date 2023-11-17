% function:     MRF_simulation
% purpose:      simulate MRF fingerprints (unaccelerated: no GPU or parallel computing in here)
%
% The theory is based on: Hancu I, Van der Maarel J, Boada F. A Model for the Dynamics of Spins 3/2 in Biological Media:Signal Loss during Radiofrequency Excitation in Triple-Quantum-Filtered Sodium MRI. Journal of Magnetic Resonance 2000; 147(2):179-191.
%
% ISTO simulations use symmetric and assymetric Tlm: Tlm(s) and Tlm(a) as described by Hancu et al.
%
% The simulations are based on the work by Arthur Magill: https://github.com/arthurmagill/isto_sim
%
% Fabian Kratzer - German Cancer Research Center (DKFZ)
% 01.03.2021 - f.kratzer@dkfz.de
%
% Haneefa's edit for 3/2 spin bSSFP simulation
% 01.09.2023 - hb3218@ic.ac.uk
clear
clc
%%
addpath('helper_functions')
%% set sample using time constant values
%{
T1=1;
T2l= 100e-3;
T2s= 5e-3;
[j0 j1 j2] = fun_calc_Js_iso(T1,T2l,T2s);
%}


%%

%prepare sample 
fluid            = struct;     
fluid.J0         = 10;          % spectral density [Hz]
fluid.J1         = 11; 
fluid.J2         = 10;
fluid.omegaQHz   = 0;           % residual quadrupolar interaction [Hz]
fluid.deltaHz    = 0;           % local field off resonance
cartilage            = struct;     
cartilage.J0         = 1225;    
cartilage.J1         = 25; 
cartilage.J2         = 24;
cartilage.omegaQHz   = 0;       
cartilage.deltaHz    = 0;  
saline             = struct;
saline.J0          = 8.9;  % spectral density [Hz]
saline.J1          = 8.9; 
saline.J2          = 8.9; 
saline.omegaQHz    = 0;    % residual quadrupolar interaction 
saline.deltaHz     = 0;

%[T1l, T1s,T2l,T2s] = fun_calc_rel_times_iso_biT1(cartilage.J0, cartilage.J1, cartilage.J2);
%%
% example data (calculate (J0,J1,J2) = (180 40 20)Hz and (220 30 15)Hz for off-resonances: -150 to 150 Hz)
w=150;
woff = cat(1,(-w:1:w).')*2*pi;           %off-resonance [Hz]*2*pi
% change the sample hete
J0 = cat(1,ones([length(woff) 1]))*saline.J0  ;         %[Hz]
J1 = cat(1,ones([length(woff) 1]))*saline.J1  ;           %[Hz]
J2 = cat(1,ones([length(woff) 1]))*saline.J2  ;             %[Hz]
wq = cat(1,ones([length(woff) 1])*0);               %[Hz]

% cat parameters
LUT_Parameter = cat(2,J0,J1,J2,wq, woff);
num_of_spins   = 1;                                       %number of isochromats
%% set FA and TE pattern
pulse_duration = 1e-3;          % pulse duration [s]
pulse_phase = zeros([1e3 1]);   % pulse phase [rad]
TE = 5.5e-3; %0.55;                      % minimal echso time [s]
B1 = 1.0;                       % relative B1 amplitude
% TE pattern
% TE_Var corresponds to TE - pulse_duration/2
%load('helper_functions/TE_Var_Na23')  
%TE_scal = 20;
%TE_Var = []
%TE_Var = (ceil((TE_Var*TE_scal)*100)/100+TE- pulse_duration/2*1e3)*1e-3; % time between pulse and acquisition: TE_scal*TE_Var + TE- d_pulse/2:
TE_Var = ones(1,1000)*(TE- pulse_duration/2); % time between pulse and acquisition: TE_scal*TE_Var + TE- d_pulse/2:
%% FA pattern
%load('helper_functions/FA_Pattern');
FA = 90;
FA_Pattern = ones(1,1000)*FA*(pi/180)*B1; %scale FA pattern
%FA_Pattern = FA_Pattern*FA_scal/100/180*pi*B1; %scale FA pattern
TR_fil = 10e-3;%13.79*1e-3; % time between sampling of k_space center and next pulse [s]
signal_length=length(TE_Var);
%% initialize  
% precalculate phase spoiling
%[cos_spoil, sin_spoil] = fun_precalc_spoiling(num_of_spins);
% signal
num_entries = size(LUT_Parameter,1);
signal = complex(zeros([signal_length num_entries]));
%% run simulation
fprintf('simulation progress: %i%% ',0)
for i_num = 1:num_entries
    % init
    Tlm = zeros([16 num_of_spins]);
    Tlm(1,:) = 1;
    Tlm(2,:) = 1;
    params.J0 = LUT_Parameter(i_num,1);
    params.J1 = LUT_Parameter(i_num,2);
    params.J2 = LUT_Parameter(i_num,3);
    params.wq = LUT_Parameter(i_num,4);
    params.woff = LUT_Parameter(i_num,5);
   
    [S,Md] = fun_calc_precessmat(params);    % precalculate precession matrices
   
    cycle_ind = 1;
    % loop over all cycles
    for FA_ind =1:length(FA_Pparametersattern)
        % excitation
        [Tlm, Tlm_history_pulse,t_pulse] = fun_rect_pulse((Tlm),params, FA_Pattern(FA_ind),pulse_phase(FA_ind),pulse_duration,1,0);
     
       
            % TE
            Tlm = fun_precess_inputmat(Tlm,TE_Var(FA_ind),S,Md);
            % ADC
            [signal_p, Tlm] = fun_acquire_point(Tlm);
            signal(cycle_ind,i_num) = signal_p;
            cycle_ind = cycle_ind+1;
     
        % TR_fil
        Tlm = fun_precess_inputmat(Tlm,TR_fil,S,Md);
        % spoiling
        %Tlm = fun_phase_spoil(Tlm,cos_spoil,sin_spoil);
       
    end
   
% show simulation status
if sum(floor(num_entries/10)*(1:9)-i_num == 0) == 1
        fprintf('%i%% ',10*i_num/floor(num_entries/10))
end
       
end
fprintf('%i%% \n',100)
%clear temp variables
%clearvars -except LUT_Parameter signal
%% show the data
off=-w:w;
figure (1)
plot(off,abs(signal(1000,:)),'LineWidth', 2); title(strcat('bSSFP for saline, FlipAngle =',num2str(FA)))
%xlim([-150,150]);
%ylim([0, 1.2*max(abs(signal(:)))]);grid on; xlabel('Off-resonance')

figure (2)
plot(abs(signal(:,100)),'LineWidth', 2); title('Na bSSFP signal evolution offer all TRs')



