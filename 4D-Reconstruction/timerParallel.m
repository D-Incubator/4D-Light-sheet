function timerParallel(outputDir,t_period_all,t_align_all,t_all,t_period,t_relativeshift,t_relativeshift_all,t_absoluteshift_all,t_finalize)

    filename = [outputDir '_timer' '.xlsx'];
    T1 = table(t_period_all,t_align_all,t_all);
    T2 = table(t_period');
    T2.Properties.VariableNames = {'t_period'};
    T3 = table(t_relativeshift');
    T3.Properties.VariableNames = {'t_relativeshift'};
%     T4 = table(t_relativeshift_all,t_absoluteshift_all,t_finalize);
%     T4.Properties.VariableNames = {'t_align'};
    writetable(T1,filename,'Sheet','all');
    writetable(T2,filename,'Sheet','get period');
    writetable(T3,filename,'Sheet','get relative shift');
%     writetable(T4,filename,'Sheet','get alignment');
    