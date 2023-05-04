function [PI_larva,v_x,v_t,spd,PI_tp,PI_tp_x,PI_tp_t,PI_tp_grad,PI_tp_grad_diff]=cal_PI_larva(x,y,t,varargin)
%this function will compute the PI for each animal at each position
v_x={};v_t={};PI_tp={};spd={};
PI_tp_x={};PI_tp_t={};PI_tp_grad={};
grad={};grad_diff={};

for i=1:2:length(varargin)
    if strcmp(varargin{i},'grad')
        grad=varargin{i+1};
    elseif strcmp(varargin{i},'grad_diff')
        grad_diff=varargin{i+1};
    end 
end 
for i=1:length(x)
    x1=x{i,1}(1:end-1,1);
    y1=y{i,1}(1:end-1,1);
    t1=t{i,1}(1:end-1,1);
    
    
    x_diff=diff(x{i,1});
    y_diff=diff(y{i,1});
    t_diff=diff(t{i,1});

    t_PI=t{i,1}(1:end-1);
    idx0=find(t_diff==0);
    idx1=find(x_diff==0&y_diff==0);
    idx3=[idx0;idx1];
    idx2=unique(idx3);

    x_diff(idx2)=[];
    y_diff(idx2)=[];
    t_diff(idx2)=[];
    t_PI(idx2)=[];
    x1(idx2)=[];
    y1(idx2)=[];
    t1(idx2)=[];
    
    if ~isempty(grad)
        grad1=grad{i,1}(1:end-1,1);
        grad1(idx2)=[];
        PI_tp_grad{i,1}=grad1;
    end 
    
    if ~isempty(grad_diff)
        grad_diff1=grad_diff{i,1}(1:end-1,1);
        grad_diff1(idx2)=[];
        PI_tp_grad_diff{i,1}=grad_diff1;
    end 

    PI_tp_x{i,1}=x1;
    PI_tp_t{i,1}=t1;
    
    vx=x_diff./t_diff;
    sp=sqrt(x_diff.^2+y_diff.^2)./t_diff;
    PI=vx./sp;

    spd{i,1}=sp;
    v_x{i,1}=vx;
    v_t{i,1}=t_PI;
    PI_tp{i,1}=PI;
    PI_larva(i,1)=mean(PI);
    %     PI_sample_size=length(dat_grouped.PI_larva);
    clearvars -except i PI_larva x y t PI_tp* v_x v_t spd grad grad_diff
end

end