
clear; close all;

load ../../../testoutput_2/Test02aMonolayerGrowth/results_from_time_0/tissuewidth.dat

plot(tissuewidth(:,1),tissuewidth(:,4), 'k')

xlabel("Time(h)")
ylabel("Tissue Width (\mu m)")

SaveAsPngEpsAndFig(-1,['Figs/ChasteMonolayerGrowth'], 12, 7/5, 12);

