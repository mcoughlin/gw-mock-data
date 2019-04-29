function[z,p,k] =  ASC_Control_Signal_Modeling(inp_filename, out_filename)
% Code to get required magnirtude response zpk model for the measured ASC Control signal
% Ref: https://git.ligo.org/NoiseCancellation/MockData/blob/master/mockdata/ASC_Models/README.md
% Prior to generating an approximate lower model, the data  smoothening and
% sharp resonances are removed. 
%
% Arguments
%   inp_filename:  Timeseries of measured control signal in .csv form [ GPS_time Counts ]
%   out_filename:  Output filename used to save .mat, .npz & .png files
%
% Usage: 
%  [z,p,k] = ASC_Control_Signal_Modeling('./Data/L1-ASC-DHARD_Y_OUT_DQ_11845440181800.csv','L1_ASC_DHARD_YAW_Model_PeakSmoothed');




%% Check & load Data
if ~nargin == 0
    x = csvread(sprintf('./Data/%s',inp_filename));
    filename = out_filename; % output file name
else
    disp('Data folder missing')
end

%% Save results
save_result = 1;


%% Bode Plot options
 BO = bodeoptions;
 BO.MagUnits='abs';
 BO.MagScale = 'log';
 BO.PhaseWrapping = 'on';
 BO.PhaseVisible = 'off';
 BO.FreqUnits = 'Hz';
 BO.Title.FontSize = 20;
 BO.XLabel.FontSize = 20;
 BO.YLabel.FontSize = 20;
 BO.TickLabel.FontSize = 20;
 BO.XLim = [1e-02 450];
 BO.YLim = [1e-3 1500];
 
%%  Sampling Rate

fs = 512;
nfft = fs*32;

 
 %% Actual & Smoothed Response
 
 [Pww, F] = pwelch(x(:,2),hanning(nfft),nfft/2,nfft,fs);
 
 Resp = sqrt(Pww);
 % Outlier removal using Hampel identifier
 Resp_hampel = hampel(Resp,10);
 % Cubic smoothing spline
 Resp_csaps = csaps(F,Resp_hampel,.9,F);
 % Savitzky-Golay filtering
 Resp_sgolay = sgolayfilt(Resp_csaps,5,9);
 % Moving average smoothing
 Resp_smooth = smooth(Resp_sgolay,50);
 

 


 

%% Get lower order zpk model for the smoothed TF response


sys = frd(Resp_smooth,F);
sys_model = fitmagfrd(abs(sys),12);


[z,p,k] = zpkdata(sys_model,'v');

disp(z)

[MAG,~] = bode(sys_model,F);
MAG = squeeze(MAG);

%% Visulalize  results
figure(33)
clf
loglog(F, Resp,'linewidth',2); axis tight; grid on
hold on 
loglog(F, Resp_smooth,'linewidth',2); axis tight; grid on
loglog(F, MAG,'linewidth',2); axis tight; grid on
xlabel('Frequency')
ylabel('Magnitude')
legend('Measured','Smoothed','Fit')

%% Save results

if save_result == 1
    
    if ~exist('Results','dir') == 1
        disp('Creating Results folder...')
    end
    
    disp('Saving Results')
    
    % Save Results
    title(filename,'interpreter','none');
    save(sprintf('./Results/Mat_models/%s.mat',filename),'z','p','k');
    saveas(gca,sprintf('./Results/Plots/%s.png',filename));
    
    % Save as .npz
    matFile = sprintf('./Results/Mat_models/%s_PeakSmoothed.mat',filename);
    pyFile =  sprintf('./Results/Py_models/%s_PeakSmoothed.npz',filename);
    system(sprintf('python convert_zpk_matlab_python.py -m %s -p %s', matFile, pyFile));
       
    
end
