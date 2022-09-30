%%% This script is for the purpose of providing outputs of empirical fragilities in an Excel file 
%%% Written By: Hossein Ebrahimian, University of Naples Federico II (UNINA)
%%% Laters Updated: 09/2022

%% Main Function

function RF2ExcelFile(output,filename,x,RF,vecDS,model,do_write_excelFile)

if do_write_excelFile == 1

data = cell(length(x)+2,7*length(vecDS));    

for j=1:length(vecDS)    

k = 7*(j-1)+1;    
data{1,k}   = vecDS{1,j};    
data{2,k}   = 'flow depth [m]';   
data{2,k+1} = 'mean-1sigma fragility';    
data{2,k+2} = 'mean fragility';    
data{2,k+3} = 'mean+1sigma fragility';  
for i=1:length(x)
    data{i+2,k}   = x(i);
    if RF.RF16(i,j)<0.0
        RF.RF16(i,j)=0.0;
    end
    if RF.RF84(i,j)>1.0
        RF.RF84(i,j)=1.0;
    end
    data{i+2,k+1} = RF.RF16(i,j); 
    data{i+2,k+2} = RF.fragility(i,j); 
    data{i+2,k+3} = RF.RF84(i,j); 
end

data{1,k+4} = 'parameters of the mean fragility';    
data{2,k+4} = 'median'; 
data{2,k+5} = 'logarithmic standard deviation';
data{2,k+6} = 'epistemic uncertainty'; 
data{3,k+4} = RF.etaIMc(j); 
data{3,k+5} = RF.betaIMc(j);
data{3,k+6} = RF.betaUF(j);

end

xlswrite([output,'\',filename,'_',model,'.xlsx'],data);

end

%% END