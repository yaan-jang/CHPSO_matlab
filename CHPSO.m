clc;clear all;close all
% % ===================================================================== %
% % Convergent Heterogenous Particle Swarm Optimization                   %
% % ===================================================================== %
% % Summary:                                                              %
% % ===================================================================== %
% %       Author: Ngaam J. Cheung                                         %
% %        Email: nyaam.ch@gmail.com                                      %
% %      Release: 1.2                                                     %
% % Release Date: June 9, 2013.                                           %
% % ===================================================================== %

% % ===================================================================== %
% % Please cite: N. J. Cheung, X.-M. Ding, H.-B. Shen, 
% %                 OptiFel: A Convergent Heterogeneous Particle Swarm 
% %                 Optimization for Takagi-Sugeno Fuzzy Modeling. 
% %                 IEEE Trans. on Fuzzy System, 
% %                 doi: 10.1109/TFUZZ.2013.2278972
% % ===================================================================== %

help CHPSO.m



NGen = 1e3;
gap = 1:(NGen/50):NGen;

% Parameters of PSO
sNP = 8;
acc = [1.49445, 1.49445];
iwt = [0.9, 0.35];  % interia weight

for funcNum = 1
    [Lb, Ub,dim] = funcRange(funcNum);  % for objFunc: twenty-three benchmark functions

    param = struct( 'ParicleNumber',           sNP,...
                        'AccelerationConstants',   acc,...
                        'InteriaWeight',           iwt,...
                        'Dimension',               dim,...
                        'Ub',                    Ub,...
                        'Lb',                    Lb);


    for run = 1:3
        rand('state',sum(100*clock))
        
        BestChart = CHPSO_func('objFunc',NGen,param,funcNum);
        eval(['CHPSO_F',num2str(funcNum),'(:,run) = BestChart;']);
        
        figure(1)
        semilogy(gap,BestChart(gap),'-ok','linewidth',1.4,'Markersize',3);
        title({strcat('\fontsize{12}\bf Function: F',num2str(funcNum));strcat('Currunt run time: ',num2str(run))});
        xlabel('\fontsize{12}\bf Generation');
        ylabel('\fontsize{12}\bf Best-so-far');
        legend('\fontsize{10}\bf CHPSO',1);
    end
    
end




