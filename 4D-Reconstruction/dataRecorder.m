function dataRecorder(outputDir,t_p,t_p_candidate,t)
    filename = [outputDir '_data' '.xlsx'];
    Final_Period = t_p;
    Potential_Period = t_p_candidate;
    Registration_Result = t;
    T1 = table(Final_Period);
    T2 = table(Potential_Period');
    T2.Properties.VariableNames = {'Potential_Period'};
    T3 = table(Registration_Result);
    T3.Properties.VariableNames = {'Registration_Result'};
    writetable(T1,filename,'Sheet','Final_Period');
    writetable(T2,filename,'Sheet','Potential_Period');
    writetable(T3,filename,'Sheet','Registration_Result');
    