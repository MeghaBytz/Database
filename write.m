clear
close
eqns=determineVdot();
edit dndt_MJC.m
FID=fopen('dndt_MJC.m','a');

% 
for j = 1:length(eqns)
    b = cellstr(eqns{j,1});
    characterEqn = char(b);
    fprintf(FID,'myV(%d)=',j);
    fprintf(FID, '%s', characterEqn);
    fprintf(FID,'%c\n',';');
end

% for j = 1:length(eqns)
% fprintf(FID, '%\n\r', eqns{j,1});
% end


fclose(FID);