function [lim]=mlt_subplt(data,series,fig_num,row,col,plt_num,fig_name,xtitle,ytitle,color,name,plot_type,varargin)
%this function will plot one subplot for each condition
%data should be a cell array that contains data for each conditon
%the first column of the data{1,1} should be the y value and the second
%column will be the std or SEM or 95% CI
%series is the x, data is y of the plot
%name is the name of each condition

% if we plot line/bar/binned bar without series, we will plot all condition in one single subplot in
% row,col,plt_num

%if we plot multiple line/bar plot with series/scatter/heatmap,
%we need to plot one subplot for each condition, we can plot them right
%after getting the data

% we plot polarplot after gathering data for all conditions because we need
% to enter the min and max values

figure(fig_num)
set(gcf,'Name',fig_name);
%% varargin
%1) min: minimum value for the data
%2) max: maximum value for the data
%3) ebar: whether we want to plot the errorbar or not:1/0
%4) legends: for naming different conditions etc
%5) ii_num: total number of conditions
%6) alpha: transparency of the scatter plot
%7) x/ytick,x/yticklabels: for heatmap
%8) title_lab:name of subplot
%9) bins is for polar plot
%10) series2: is for imagesc(X,Y,C) for the heatmap
mi=0;ma=10;std=0;legends={};ii_num=0;alpha=0.5;title_lab='';series2=[];
line_type='-'; marker_size=0.3;
rlims=[]; thetalims=[];
y_lim=[];x_lim=[];
ydir=0; lw=1;
bt_CI=0;% plot the 95% CI from bootstrap or not
bar_color=[];
yline1=[];
labels={};
for i=1:2:length(varargin)
    if strcmp(varargin{i},"min")
        mi=varargin{i+1};
    elseif strcmp(varargin{i},"max")
        ma=varargin{i+1};
    elseif strcmp (varargin{i},"ebar")
        std=varargin{i+1};
    elseif strcmp (varargin{i},"legends")
        legends=varargin{i+1};
    elseif strcmp (varargin{i},"ii_num")
        ii_num=varargin{i+1};
    elseif strcmp (varargin{i},"alpha")
        alpha=varargin{i+1};
    elseif strcmp (varargin{i},"xtick")
        xt=varargin{i+1};
    elseif strcmp (varargin{i},"ytick")
        yt=varargin{i+1};
    elseif strcmp (varargin{i},"xticklabel")
        xtl=varargin{i+1};
    elseif strcmp (varargin{i},"yticklabel")
        ytl=varargin{i+1};
    elseif strcmp (varargin{i},"title")
        title_lab=varargin{i+1};
    elseif strcmp (varargin{i},"bins")
        bins=varargin{i+1};
    elseif strcmp (varargin{i},"series2")
        series2=varargin{i+1};
    elseif strcmp (varargin{i},"line_type")
        line_type=varargin{i+1};
    elseif strcmp (varargin{i},"marker_size")
        marker_size=varargin{i+1};
    elseif strcmp(varargin{i},'rlim')
        rlims=varargin{i+1};
    elseif strcmp(varargin{i},'tlim')
        thetalims=varargin{i+1};
    elseif strcmp(varargin{i},'ylim')
        y_lim=varargin{i+1};
    elseif strcmp(varargin{i},'xlim')
        x_lim=varargin{i+1};
    elseif strcmp(varargin{i},'ydir')
        ydir=varargin{i+1};
    elseif strcmp(varargin{i},'lw') %linewidth
        lw=varargin{i+1};
    elseif strcmp(varargin{i},'bt_CI') %linewidth
        bt_CI=varargin{i+1};
    elseif strcmp(varargin{i},'bar_color') %linewidth
        bar_color=varargin{i+1};
    elseif strcmp(varargin{i},'yline') %linewidth
        yline1=varargin{i+1};
    elseif strcmp(varargin{i},'label') %linewidth
        labels=varargin{i+1};
    end
end
%% plot line plot with or without patch
if strcmp(plot_type,"line")
    subplot(row,col,plt_num)
    hold on
    %eliminate the nan value from the array
    len=[];idx={};
    for i=1:ii_num

        idx{i,1}=find(isnan(data{i,1}(:,1))|isinf(data{i,1}(:,1)))
        len(i,1)=length(idx{i,1})
    end
    [M,I]=max(len);
    if ~isempty(idx)
        series(idx{I,1},:)=[];
    end
    %if we want to just simply plot line without any std or standard error
    for i=1:ii_num

        if ~isempty(idx)
            data{i,1}(idx{I,1},:)=[];
        end
        if std==1

            patch([series;flipud(series)],[data{i,1}(:,1)-data{i,1}(:,2);flipud(data{i,1}(:,1)+data{i,1}(:,2))],...
                color(i,:),'EdgeColor',color(i,:),'FaceAlpha',0.4,'EdgeAlpha',0.4);
        elseif bt_CI==1
            %plot the 95% CI from bootstrap data-->it can be asymetric
            patch([series;flipud(series)],[data{i,1}(:,1)+data{i,1}(:,3);flipud(data{i,1}(:,1)+data{i,1}(:,2))],...
                color(i,:),'EdgeColor',color(i,:),'FaceAlpha',0.4,'EdgeAlpha',0.4);
        end
        if ~isempty(color)
            p(i)= plot(series,data{i,1}(:,1),'-','Color',color(i,:),'LineWidth',1);
        else
            p(i)= plot(series,data{i,1}(:,1),'-','LineWidth',1);
        end
        hold on
    end
    ax=gca;
    axis tight
    title(title_lab,'FontSize',28);
    if ~isempty(name)
        legend(p,name,'FontSize',20);
    end
    %flip the y axis
    if ydir==1
        axis ij
    end
    if ~isempty(y_lim)
        ylim(y_lim);
    end
    if ~isempty(x_lim)
        xlim(x_lim);
    end

    ax.FontSize=20;
    ax.FontWeight='bold';
    xlabel(xtitle,'FontSize',24);
    ylabel(ytitle,'FontSize',24);

    hold off

elseif strcmp(plot_type,"multiple lines")
    %this is for plotting multiple lines for the same condition in one
    %subplot, for example, across the time, the number of larvae in three
    %different regions

    subplot(row,col,plt_num)
    hold on
    idx=find(isnan(data(:,1))|isinf(data(:,1)));
    data(idx,:)=[];
    series(idx)=[];
    clear idx
    for j=1:width(data)
        if line_type=='-'
            p(j,1)=plot(series,data(:,j),'-','LineWidth',lw,'Color',color(j,:));
        elseif line_type=='.'
            p(j,1)=plot(series,data(:,j),'.','Color',color(j,:),'MarkerSize',0.5);
        end
    end
    title(name,'FontSize',28)
    ax=gca;
    ax.FontSize=20;
    ax.FontWeight='bold';
    xlabel(xtitle,'FontSize',24);
    ylabel(ytitle,'FontSize',24);
    if ~isempty(legends)
        legend(p,legends,'FontSize',20); %legends need to be a cell array
    end

    if ~isempty(y_lim)
        ylim(y_lim);
    end
    axis tight
    hold off

    if plt_num==ii_num
        for i=1:ii_num
            figure(fig_num);
            l(i,1)=subplot(row,col,i);
        end
        linkaxes(l,'y');
    end

    %% plot bar plot
elseif strcmp (plot_type,"bar")
    %depends on the data type there are two different plots
    %1)the input data is a matrix similar to the output of the mean value
    %of preferential index
    if isempty(series)
        if iscell(data)
            %2) if this data contains onlycells, one for the mean value one for the std
            %each have the datapoints for all conditions
            subplot(row,col,plt_num);

            hBar = bar(data{1,1},'grouped','EdgeColor','none');
            [ngroups,nbars] = size(data{1,1});
            x = nan(nbars,ngroups);
            for i = 1:nbars
                x(i,:) = hBar(i).XEndPoints;
            end
            % Plot the errorbars
            hold on
            if std==1
                errorbar(x',data{1,1},data{1,2},'k','linestyle','none');
            elseif bt_CI==1
                errorbar(x',data{1,1},data{1,3},data{1,2},'k','linestyle','none');
            end
            if ~isempty(bar_color)
                for k=1:height(bar_color)
                    hBar(k).FaceColor=bar_color(k,:);
                end
            end

            legend(legends,'FontSize',20);
            l=height(name);
            xlim([0 l+1]);
            xticks(0:1:l+1);
            name=[' ';name];
            xticklabels(name);
            ax=gca;
            %             ax.FontSize=20;


            if ~isempty(y_lim)
                ylim(y_lim);
            end

            if ~isempty(yline1)
                for len1=1:length(yline1)
                    yline(yline1(len1),'r--','LineWidth',1.5);
                end
            end

            axis square;
            xlabel(xtitle,'FontSize',24);
            ylabel(ytitle,'FontSize',24);
            title_lab
            title(title_lab,'FontSize',28);
            hold off;

        elseif ismatrix(data)


            subplot(row,col,plt_num)
            title(title_lab,'FontSize',28);
            hold on
            hBar = bar(data(:,1),'EdgeColor','none');
            if std==1
                errorbar(data(:,1),data(:,2),'.k')
            elseif bt_CI==1
                xpos=[1:length(data(:,1))]';
                errorbar(xpos,data(:,1),data(:,3),data(:,2),'.k')
            end

            xlim([0 ii_num+1]);
            xticks(0:1:ii_num+1);
            name=[' '; name];
            xticklabels(name);
            ax=gca;
            ax.FontSize=20;
            ax.FontWeight='bold';

            xlabel(xtitle,'FontSize',24);
            ylabel(ytitle,'FontSize',24);
            hBar.FaceColor='flat';
            for i=1:ii_num
                hBar.CData(i,:)=color(i,:);
            end

            if ~isempty(y_lim)
                ylim(y_lim);
            end
            axis square;
            hold off
        end
    elseif ~isempty(series)
        %for other bar data, they will contain the same number of
        %conditions as the name does, each cell array will have two
        %columns, one for mean value or count etc, one for std
        subplot(row,col,plt_num)
        hold on
        hBar=bar(series,data(:,1),'FaceColor',color(plt_num,:),'EdgeColor','none');
        if std==1
            errorbar(series,data(:,1),data(:,2),'.k')
        elseif bt_CI==1
            errorbar(series,data(:,1),data(:,3),data(:,2),'.k')
        end
        title(name,'FontSize',28);
        ax=gca;
        ax.FontSize=20;
        ax.FontWeight='bold';

        xlabel(xtitle,'FontSize',24);
        ylabel(ytitle,'FontSize',24);
        if ~isempty(y_lim)
            ylim(y_lim);
        end
        if ~isempty(x_lim)
            xlim(x_lim);
        end
        hold off

        if plt_num==ii_num
            for i=1:ii_num
                figure(fig_num);
                l(i,1)=subplot(row,col,i);
            end
            linkaxes(l,'y');
        end

    end
elseif strcmp(plot_type,'overlapped bar')
    %% this will plot the data vs series from two groups: exp vs ctr
    %the input data should be a cell contains the data you want to plot in
    %the same histogram
    % for each cell in hte data cell array, it can have two columns, one
    % for mean and one for sem
    subplot(row,col,plt_num)
    hold on
    for i =1:length(data)
        bar(series,data{i,1}(:,1),'FaceColor',color(i,:),'FaceAlpha',0.7,'EdgeColor','none');
        if std==1
            errorbar(series,data{i,1}(:,1),data{i,1}(:,2),'.k')
        elseif bt_CI==1
            errorbar(series,data{i,1}(:,1),data{i,1}(:,3),data{i,1}(:,2),'.k')
        end
    end

    title(name);
    xlabel(xtitle);
    ylabel(ytitle);
    ax=gca;
    ax.FontSize=20;
    ax.FontWeight='bold';
    if ~isempty(legends)
        legend(legends,'FontSize',20);
    end
    if plt_num==ii_num
        for i=1:ii_num
            figure(fig_num);
            l(i,1)=subplot(row,col,i);
        end
        linkaxes(l,'y');
    end
    if ~isempty(y_lim)
        ylim(y_lim);
    end

    if ~isempty(x_lim)
        xlim(x_lim);
    end
    axis square
    hold off
elseif strcmp(plot_type,'stacked bar')
    %% plot stacked bar giving the data and series
    subplot(row,col,plt_num)
    hold on
    axis square
    b=bar(series,data,'stacked','EdgeColor','none');
    for k=1:width(data)
        b(k).CData=color(k,:);
    end 
    %     title(name,'FontSize',28);
    ax=gca;
    ax.FontSize=24;
    ax.FontWeight='bold';
    xlabel(xtitle,'FontSize',24);
    ylabel(ytitle,'FontSize',24);
    title(title_lab);
   

    if ~isempty(legends)
        legend(legends);
    end

    if ~isempty(y_lim)
        ylim(y_lim);
    end

    xlim([0 ii_num+1]);
    xticks(0:1:ii_num+1);
    name=[' '; name];
    xticklabels(name);
    % label the data inside the bar 
    if ~isempty(labels)
        xt=[1:ii_num]';
        barbase=cumsum([zeros(size(data,1),1) data(:,1:end-1)],2);
        labelpos=data/2+barbase;

        for k1=1:size(data,1)

            text(xt(k1)*ones(1,size(data,2)),labelpos(k1,:),labels(2*k1-1:2*k1),'HorizontalAlignment','center','FontSize',20,'FontWeight','Bold');
        end 
    end 
    
    hold off
elseif strcmp(plot_type,"fitted line")

    subplot(row,col,plt_num)
    hold on
    tbl=table(data,series);
    mdl=fitlm(tbl,'linear');
    plot(mdl,'Marker','.','Markersize',3,'Color',color(plt_num,:));
    title(name,'FontSize',28);
    ax=gca;
    ax.FontSize=20;
    ax.FontWeight='bold';

    xlabel(xtitle,'FontSize',24);
    ylabel(ytitle,'FontSize',24);
    axis tight
    axis square
    hold off

elseif strcmp (plot_type,"polar plot")
    clear series
    series=deg2rad(-180:bins:180);
    series_sm=deg2rad(-180+bins/2:bins:180-bins/2);


    for i=1:length(data)
        subplot(row,col,i,polaraxes)

        arry=data{i,1}(:,1)-mi;
        polarhistogram('BinEdges',series,'BinCounts',arry,'FaceColor',color(i,:),'EdgeColor',[1 1 1]);
        hold on



        if std==1
            upper=data{i,1}(:,1)+data{i,1}(:,2);
            lower=data{i,1}(:,1)-data{i,1}(:,2);
        elseif bt_CI==1
            upper=data{i,1}(:,1)+data{i,1}(:,2);
            lower=data{i,1}(:,1)+data{i,1}(:,3);
        end

        if std==1||bt_CI==1

            for j=1:length(upper)
                if  ~isnan(upper(j))
                    polarplot([series_sm(j),series_sm(j)],[upper(j)-mi,lower(j)-mi],'k-','Linewidth',0.7);
                end
            end
        end



        if ~isempty(name)
            title(name(i),'FontSize',20)
        end

        thetaticks([0:45:360])
        thetaticklabels([0 45 90 135 180 -135 -90 -45])
        title(name(i),'FontSize',20);
        rlim([0 ma-mi]);
        rt=rticks;
        rt=rt+mi;
        rticklabels(round(rt,2));

        if ~isempty(thetalims)
            thetalim(thetalims);

        end
        hold off
    end

elseif strcmp(plot_type,"overlapped polar plot")
    clear series
    series=deg2rad(-180:bins:180);
    series_sm=deg2rad(-180+bins/2:bins:180-bins/2);
    subplot(row,col,plt_num)

    for i=1:length(data)
        length(data{i,1})
        length(series_sm)
        arry=data{i,1}(:,1)-mi;
        polarhistogram('BinEdges',series,'BinCounts',arry,'FaceColor',color(i,:),'EdgeColor',[1 1 1]);
        hold on
        if std==1
            upper=data{i,1}(:,1)+data{i,1}(:,2);
            lower=data{i,1}(:,1)-data{i,1}(:,2);
        elseif bt_CI==1
            upper=data{i,1}(:,1)+data{i,1}(:,2);
            lower=data{i,1}(:,1)+data{i,1}(:,3);
        end
        if std==1||bt_CI==1
            for j=1:length(upper)
                if  ~isnan(upper(j))
                    polarplot([series_sm(j),series_sm(j)],[upper(j)-mi,lower(j)-mi],'k-','Linewidth',0.7);
                end
            end
        end

    end
    rlim([0 ma-mi]);
    rt=rticks;
    rt=rt+mi;
    rticklabels(round(rt,2));

    thetaticks([0:45:360]);
    thetaticklabels([0 45 90 135 180 -135 -90 -45]);
    if ~isempty(thetalims)
        thetalim(thetalims);
    end
    if ~isempty(name)
        title(name,'FontSize',20)
    end
    hold off
elseif strcmp(plot_type,"scatter")
    subplot(row,col,plt_num)
    hold on;
    scatter1=scatter(series,data,marker_size,color(plt_num,:));
    scatter1.MarkerEdgeAlpha=alpha;
    title(name(plt_num),'FontSize',28);
    ax=gca;
    ax.FontSize=20;
    ax.FontWeight='bold';
    xlabel(xtitle,'FontSize',24);
    ylabel(ytitle,'FontSize',24);
    if ~isempty(legends)
        legend(legends,'FontSize',20);
    end
    if ~isempty(y_lim)
        ylim(y_lim);
    end

    if ~isempty(x_lim)
        xlim(x_lim);
    end

    axis tight
    hold off
elseif strcmp(plot_type,"heatmap")
    subplot(row,col,plt_num)
    hold on;
    imagesc(series,series2,data);
    title(name,'FontSize',18);
    ax=gca;
    ax.FontSize=18;
    ax.FontWeight='bold';
    xlabel(xtitle,'FontSize',18);
    ylabel(ytitle,'FontSize',18);
    colorbar;
    lim=caxis;
    axis tight
    axis square
    hold off
elseif strcmp(plot_type,"boxplot")
    subplot(row,col,plt_num)
    hold on;
    boxplot(data,series,"PlotStyle",'compact');
    if ~isempty(yline1)
        for len1=1:length(yline1)
            yline(yline1(len1),'r--','LineWidth',1.5);
        end
    end
    if ~isempty(y_lim)
        ylim(y_lim);
    end
    axis square
    xlabel(xtitle,'FontSize',18);
    ylabel(ytitle,'FontSize',18);
    title(name,'FontSize',18);
    ax=gca;
    ax.FontSize=20;
    ax.FontWeight='bold';

    hold off

end
end