function []=add_opto_grad_choreo(tracker_num,genotype,condition,filename,opto_cond)
%this function will use the fitting curve to compute the opto gradient
%for each x pos--Only used for optogenetic exp

%opto_cond needs to be ["C13";"C14";"C15";"C16";"C17"]
%set output folder
outdir=fullfile("/project/6010970/screen/olfactory_output/choreography",tracker_num,genotype,condition)

%% Loading the coefficient for the gradient fitting curve
% we need to check the coef is for which conditions (code is in python, but data collected from the excel sheet)

fileID=fopen("/project/6010970/screen/olfactory_output/light_gradient_coef.txt",'r');
delimiter = ' ';
startRow = 0;
formatSpec = '%f%f%f%f';
coef= textscan(fileID, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
coef=cell2mat(coef)
fclose(fileID);
gradient=["C13";"C14";"C15";"C16";"C17"];

%% Load the existing data and use x postion for the opto gradient
dat_loc=fullfile(outdir,filename)
load(dat_loc);
k=find(matches(gradient,opto_cond));
coef_x=coef(k,:);
if k==5
    syms f(x)
    f(x)=(coef_x(1)/(1+exp(-coef_x(3)*(x-coef_x(2)))))+coef_x(4);
    ff=matlabFunction(f);
    g=diff(f);
    gg=matlabFunction(g);
    clear f g

    disp(append("The fitted curve for",opto_cond,": ", num2str(coef_x(1)),"/(1+e^(-",num2str(coef_x(3)),"*(x-",num2str(coef_x(2)),"))))+",num2str(coef_x(4))));
else
    syms f(x)
    f(x)=coef_x(1)*x^3+coef_x(2)*x^2+coef_x(3)*x+coef_x(4);
    ff=matlabFunction(f);
    g=diff(f);
    gg=matlabFunction(g);
    clear f g

    disp(append("The fitted curve for ",opto_cond, ": ",num2str(coef_x(1)),' x^3+',num2str(coef_x(2)),' x^2+ ',num2str(coef_x(3)),' x+',num2str(coef_x(4))));
end
%load('C:\Users\feihu\OneDrive - McGill University\matlab\Analysis based on larvae\data.mat')
for j=1:length(dat_grouped.x)
    x1=dat_grouped.x{j,1};
    grad{j,1}=ff(x1);
    grad_diff{j,1}=gg(x1);
end
%% save data to the var dat
dat_grouped.grad=grad;
dat_grouped.grad_diff=grad_diff;
%% plot the gradient and gradient diff curves just to confirm
g(:,1)=vertcat(dat_grouped.grad{:});
gd(:,1)=vertcat(dat_grouped.grad_diff{:});
x2=vertcat(dat_grouped.x{:});
mlt_subplt(g,x2,1,1,2,1,'','X position (mm)','Light Intensity (uW/mm^2)',[1,0,0],opto_cond,'multiple lines','line_type','.');
mlt_subplt(gd,x2,1,1,2,2,append("Fitted curves for light gradient and intensity for ",opto_cond),'X position (mm)','Light Intensity (uW/mm^2)',[0,0,0],opto_cond,'multiple lines','line_type','.');
%%
save_all_figures(outdir);

%% save the data file for each condition

filename1=fullfile(outdir,'data_grad.mat')
if ~isfolder(outdir)
    mkdir(outdir)
end

if isfile(filename1)
    delete(filename1);
end
save(filename1,'dat_grouped');


end