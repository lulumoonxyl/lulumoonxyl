function []=add_opto_grad_JB(tracker_num,genotype,condition,filename,opto_cond)
%this function will use the fitting curve to compute the opto gradient
%for each x pos--Only used for optogenetic exp

%opto_cond needs to be ["C13";"C14";"C15";"C16";"C17"]
%set output folder
outdir=fullfile("/project/6010970/screen/olfactory_output/JB_JAABA",tracker_num,genotype,condition)

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
dat_loc=fullfile(outdir,filename);
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

    disp(append("The fitted curve for ",opto_cond,": ", num2str(coef_x(1)),"/(1+e^(-",num2str(coef_x(3)),"*(x-",num2str(coef_x(2)),"))))+",num2str(coef_x(4))));
else
    syms f(x)
    f(x)=coef_x(1)*x^3+coef_x(2)*x^2+coef_x(3)*x+coef_x(4);
    ff=matlabFunction(f);
    g=diff(f);
    gg=matlabFunction(g);
    clear f g

    disp(append("The fitted curve for ",opto_cond,": ",num2str(coef_x(1)),' x^3+',num2str(coef_x(2)),' x^2+ ',num2str(coef_x(3)),' x+',num2str(coef_x(4))));
end
%load('C:\Users\feihu\OneDrive - McGill University\matlab\Analysis based on larvae\data.mat')
for j=1:length(dat.x)
    x1=dat.x{j,1};
    xspine=dat.xspine{j,1};
    turn_x=dat.turn_x{j,1};
    grad_x{j,1}=ff(x1);
    grad_turn_x{j,1}=ff(turn_x);
    grad_xspine{j,1}=ff(xspine);

    grad_t0{j,1}=grad_x{j,1}(dat.t0_idx{j,1});
    grad_t1{j,1}=grad_x{j,1}(dat.t1_idx{j,1});
    grad_r0{j,1}=grad_x{j,1}(dat.r0_idx{j,1});
    grad_r1{j,1}=grad_x{j,1}(dat.r1_idx{j,1});

    grad_turn_t0_head{j,1}=grad_xspine{j,1}(dat.t0_idx{j,1},1);
    grad_turn_t1_head{j,1}=grad_xspine{j,1}(dat.t1_idx{j,1},1);

    grad_diff_x{j,1}=gg(x1);
    grad_diff_turn_x{j,1}=gg(turn_x);
    grad_diff_xspine{j,1}=gg(xspine);

    grad_diff_t0{j,1}=grad_diff_x{j,1}(dat.t0_idx{j,1});
    grad_diff_t1{j,1}=grad_diff_x{j,1}(dat.t1_idx{j,1});
    grad_diff_r0{j,1}=grad_diff_x{j,1}(dat.r0_idx{j,1});
    grad_diff_r1{j,1}=grad_diff_x{j,1}(dat.r1_idx{j,1});
    
    grad_diff_turn_t0_head{j,1}=grad_diff_xspine{j,1}(dat.t0_idx{j,1},1);
    grad_diff_turn_t1_head{j,1}=grad_diff_xspine{j,1}(dat.t1_idx{j,1},1);
end
%% save data to the var dat
dat.grad_x=grad_x;
dat.grad_turn_x=grad_turn_x;
dat.grad_xspine=grad_xspine;
dat.grad_t0=grad_t0;
dat.grad_t1=grad_t1;
dat.grad_r0=grad_r0;
dat.grad_r1=grad_r1;
dat.grad_turn_t1_head=grad_turn_t1_head;
dat.grad_turn_t0_head=grad_turn_t0_head;

dat.grad_diff_x=grad_diff_x;
dat.grad_diff_turn_x=grad_diff_turn_x;
dat.grad_diff_xspine=grad_diff_xspine;
dat.grad_diff_t0=grad_diff_t0;
dat.grad_diff_t1=grad_diff_t1;
dat.grad_diff_r0=grad_diff_r0;
dat.grad_diff_r1=grad_diff_r1;
dat.grad_diff_turn_t1_head=grad_diff_turn_t1_head;
dat.grad_diff_turn_t0_head=grad_diff_turn_t0_head;
%% save the data file for each condition

filename1=fullfile(outdir,'data_grad.mat');
if ~isfolder(outdir)
    mkdir(outdir)
end

if isfile(filename1)
    delete(filename1);
end
save(filename1,'dat','-v7.3');


end