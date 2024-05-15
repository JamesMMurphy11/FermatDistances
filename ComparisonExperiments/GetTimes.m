%% Extract Times from the Results Structs

for i=1:10
    FD_Times(i)=sum(cell2mat(FD{i}.Times));
   
end

FD_Times_Avg=mean(FD_Times)
FD_Times_SD=sqrt(var(FD_Times))

for i=1:10
    ED_Times(i)=sum(cell2mat(ED_Rescaled{i}.Times));
end

ED_Times_Avg=mean(ED_Times)
ED_Times_SD=sqrt(var(ED_Times))
