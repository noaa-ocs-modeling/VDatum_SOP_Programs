% Example of using known modeled datum errors at gages 
% for estimating how these errors correlation coefficient  
% depends on separation (of any kind) between the gages, and
% for testing hypotheses about how this coefficient can be 
% approximated with some prescribed function(s)  
%----------------------------------------------------------
% See function correlation_estimator below the example
%----------------------------------------------------------
% variables:
% m, int32 - number of gages
% mdist, m x m double - matrix of gage2gage distances returned by wdist 
%              (measured on the mesh along the cell edges) or travel times
% cid, 1 x m int32 - string of gages' IDs
% mhh,mhw,mlw,mll, 1 x m single - observed datums
% lv, 8 x m double - modeled datums at gages as in TAD output, 
%           column: 1-#, 2-mhhw,3-mhw, 4-msl,5,6-dtl/mtl,7-mlw,8-mllw  
% whatelse, m x m double - matrix of whatever factor is being tested  
%     for affecting the error correlation (for instance, coef(i,j) - 
%     correlation of tidal envelopes), set to ones(m) is none used
% eestn, 4 x m double - rows of known model errors for mhh,mhw,mlw,mll
% xmax - max distance/time to look for any meaninfull correlation
% nb - number of bins to divide [0 xmax] into
% xbins, 1 x nb double - bin edges on the right side
%------ output of correlation_estimator -------------
% nnn, 1 x nb int32 - number of products in a bin
% xxx, 1 x nb double - average gage2gage distance in a bin
% yyy, 4 x nb double - error correlation coefficient estimate
%                                      in a bin, per datum
%---------------------------------------------------
% L - correlation distance/time as estimated from fitting yyy vs xxx
%%%%%%%%%%%%%%%   elena.tolkova@noaa.gov   %%%%%%%%%%%%

%*** gage order must be the same in all data structures ***

clear variables
its_shortest=1;
if its_shortest
    pathin='wdist/';
    xlabl='km';
    xmax=600; %largest distance at which correlation calculation stops
else
    pathin='wtime/';
    xlabl='min';
    xmax=500;  %largest travel time at which correlation calculation stops
end

path_grid='./';
eval(['load ' path_grid 'pmoe_datum_301'])

cid=[pmoe_datum.id];     % station id
nodes=int32([pmoe_datum.node]);  % grid node # at station location
m=length(nodes)
% -------------------- observed datums at gages -----------------
mhh= [pmoe_datum.omhhw];
mhw= [pmoe_datum.omhw];
mlw= [pmoe_datum.omlw];
mll= [pmoe_datum.omllw];
% -------------------- calculated datums at gages -----------------
load a79_rerun58_BC_TADxx
lv=testxx(:,nodes);
clear testxx
%---------------- convert model datums from Model Zero to MSL ----------
lvmsl=lv(3,:);
msk=(lvmsl<-999)|isnan(lvmsl);
for i=1:7
    lv(i,:)=lv(i,:)-lvmsl;
end
lv(1:7,msk)=NaN;
clear msk
% --------------------- known datum errors ----------------------------
eestn(1:4,1:m)=NaN;    % model datum error at stations
eestn(1,:)=mhh-lv(1,:);
eestn(2,:)=mhw-lv(2,:);
eestn(3,:)=mlw-lv(6,:);
eestn(4,:)=mll-lv(7,:);
% ------------- remove gages with extra large errors -------------
estnmx=max(abs(eestn));
msk=estnmx>0.4;
eestn(:,msk)=[];
cid(msk)=[];
nodes(msk)=[];
m=length(nodes)
% -------------------- gage 2 gage distances, wdist -----------------
mdist(1:m,1:m)=NaN;
for i=1:m
    iid=cid(i);
    infile=[pathin 'station' int2str(iid) '.nc'];
    dis=ncread(infile,'dis'); % water distance from i-th station
    mdist(i,:)=dis(nodes);
end
msk=(mdist>2000);  % dry nodes have mdist=Inf
mdist(msk)=NaN;
clear dis
% ---------------- bins for averaging corr products -------------
nb=30; 
xbins=logspace(0,1,nb); xbins=(xmax/10)*(xbins-1)+1;

whatelse=ones(m);

[nnn,xxx,yyy]=correlation_estimator(m,eestn,mdist,whatelse,xbins,1,its_shortest,xlabl);
save Alaska_time xmax nnn xxx yyy

% mark bins containing fewer than nmin entries for removal in ex,ey
clear msk
minbin=1;
msk=(nnn<minbin);
xxx(msk)=NaN;
yyy(1:4,msk)=NaN;

figure
axes('Position',[0.08 0.15 0.88 0.8]);
hold on
for jj=1:4
    plot(xxx,yyy(jj,:),'.','markersize',20)
%     ex=[ex, xxx(1:K)];
%     ey=[ey, yyy(jj,1:K)];
end
ex=xxx;
ey=sum(yyy)/4;
plot(ex,ey,'^','markersize',10,'markerfacecolor','r')

clear msk
msk=isnan(ex);
ex(msk)=[];
ey(msk)=[];

xlim([0 xmax])
grid on
gca.FontSize=12;
xlabel(xlabl,'fontsize',14);

% fit correlation coefficient estimate (ey vs ex) 
Ke=length(ex);
LL=0:0.1:100;
nL=length(LL);
clear r rr
r(1:nL)=0;
rr(1:Ke)=0;
for i=2:nL
    L=LL(i);
    for j=1:Ke
        rr(j)=abs(ey(j)-exp(-sqrt(ex(j)/L)));
%        rr(j)=abs(ey(j)-L/(ex(j)+L));
    end
    r(i)=sum(rr);
end
[tmp iL]=min(r(2:end))
L=ceil(LL(iL))
xr=0:0.1:500;
for j=1:length(xr)
    rr(j)=exp(-sqrt(xr(j)/L));
%    rr(j)=L/(xr(j)+L);
end
hold on
plot(xr,rr,'linewidth',2)
xlim([0 xmax/2.9412])
ylim([-0.1 1])
yticks([0:0.2:1])
grid on
legend('mhh','mhw','mlw','mll','mean',['e(-x/L),L=', int2str(L),xlabl])
title('error correlation estimate vs mdist')
%----------------------------------------------------------
% this function uses known modeled datum errors at gages 
% for estimating how these errors correlation coefficient  
% depends on distance (of any kind) between the gages, and
% for testing hypotheses about how this coefficient can be 
% approximated with some prescribed function(s)    

function [nnn,xxx,yyy]=correlation_estimator(m,eestn,dist,whatelse,xbins,iplotprod,its_shortest,xlabl)

    nb=length(xbins);
    if iplotprod  % plot error products for m(m-1)/2 gage pairs, mhhw
    
        figure
        ax=axes('Position',[0.04 0.18,0.85,0.75]);
        yyaxis left
        hold on
        xlim([0 xbins(end)])

        for jj=1:4
            z=eestn(jj,:);
            a=mean(z);
            estn =z-a;
            m_error = std(estn);
            estn=estn/m_error;
            ecorr=estn'*estn;
            if jj<3
                c='b';
            else
                c='g';
            end

            for i=1:m
                plot(dist(i,i:m),ecorr(i,i:m),'.','color',c,'markersize',9)
            end  
        end
        xticks(xbins)
        xlb=cell(1,nb);
        for i=6:3:nb
            xlb{i}=int2str(round(xbins(i)));
        end
        ax=gca;
        ax.XTickLabels=xlb;
        xtickangle(45)
        ax.TickLength=[0.05 0.005]
        %yticks([-20:10:30])
        ax.FontSize=12;
    end

%                      %  placeholders for:
    xxx(1:nb)=NaN; % average g2g distance within a bin
    nnn(1:nb)=int32(0); % number of values within a bin
    yyy(1:4,1:nb)=NaN; % mean ecorr within a bin, for each of 4 datums

    for jj=1:4
        z=eestn(jj,:);
        a=mean(z);
        estn = z-a;
        m_error = std(estn);
        estn=estn/m_error;
    % ecorr(i,j)=product of normalized jj-th datum errors at i-th and j-th gages
        ecorr=estn'*estn; 
    % divide ecorr element-wise by hypothesized correlation factor whatelse
        ecorr=ecorr./whatelse;        
    % remove duplicates and reshape into a single row
        for i=2:m
            dist(i,1:(i-1))=NaN;
            ecorr(i,1:(i-1))=NaN;
        end
        xx=reshape(dist,1,m*m);
        yy=reshape(ecorr,1,m*m);
        msk=isnan(xx);
        xx(msk)=[];
        yy(msk)=[];
    %                      %  placeholders for:
        xxx(1:nb)=NaN; % average g2g distance within a bin
        nnn(1:nb)=int32(0); % number of values within a bin
        yyy(jj,1:nb)=NaN; % sum of ecorr within a bin, jj-th datum
        for i=1:nb
            clear msk
            msk=xx<xbins(i);
            k=sum(msk);
            if k>0
                nnn(i)=k;
                xxx(i)=sum(xx(msk))/k;
                yyy(jj,i)=sum(yy(msk))/k;
                xx(msk)=[];
                yy(msk)=[];
            end
        end
    end
    
    if iplotprod       
        yyaxis right
        xlim([0 xbins(end)])
        bar(xxx,nnn,'edgeColor','none');
        xlabel(xlabl,'fontsize',14)
    end

end
