% toylandscape.m
% --------------
%{
construct a toy model to illustrate the method

%}
clear 
% close all
clc

% w=0.5;
% a=1;
% b=2;
% c=0;
% f=@(x) w*x + a*cos(x) + b*cos(2*x) + c*cos(3*x);

a=2;%-4;%
% b=-5;
% c=1;
% d=0;

b=3;b_max=4;
b_=0;

c=-4;
c_=0.1;
d=1;


% -- limits
xl = [-2.5 2.5];
yl = [-20 15];

% -- varying parameters
Npar=200;
t = linspace(-1,1,Npar);

A = linspace(a,-a,Npar);
AF = @(x) -a*x;
B = [linspace(b,b_max,Npar/2), linspace(b_max,b,Npar/2)];
BF = @(x) (b_max-b)*cos(pi*x/range(t)) + b;% input from -1 to 1
C_ = linspace(c_/1.5,-c_,Npar);
CF_ = @(x) -c_*5/6*x - 1/6*c_; 
S = [linspace(1,2,Npar/2), linspace(2,1,Npar/2)];% scaling factor
S = S.^2;
SF = @(x) cos(pi*x/range(t)) +1;

% -- plot

% -- x values
Nx = 200;
X = linspace(xl(1),xl(2),Nx);

f=@(x) a.*x + b.*(x-b_).^2 + c.*(x-c_).^4 + d.*x.^6;
figure
plot(X,f(X))
ylim(yl)
xlim(xl)

% return
% --animation
figure
% hold on
% plot(X,f(X))
ylim(yl)
xlim(xl)



for n=1:Npar
%     f=@(x) a.*x + b.*x.^2 + c.*x.^4 + d.*x.^6;
    f=@(x) S(n).*(A(n).*x + B(n).*x.^2 + c.*(x-C_(n)).^4 + d.*x.^6);
    f_=@(x) SF(t(n)).*(AF(t(n)).*x + BF(t(n)).*x.^2 + c.*(x-CF_(t(n))).^4 + d.*x.^6);%c.*(x-CF_(t(n))).^4
    plot(X',[f(X)',f_(X)'])
    ylim(yl)
    xlim(xl)
    drawnow
    pause(0.01)
end

% xltrim =pi/5;
% xl = [-2*pi+xltrim,2*pi-xltrim];

% [NN,XX]=meshgrid(1:Npar,X);
% 
% F = @(t,x) S(t).*(A(t).*x + B(t).*x.^2 + c*(x-C_(t)).^4 + d*x.^6);
% DF = @(t,x) S(t).*(A(t) + 2*B(t).*x + 4*c*(x-C_(t)).^3 + 6*d*x.^5);
% DDF = @(t,x) S(t).*(2*B(t) + 12*c*(x-C_(t)).^2 + 30*d*x.^4);
% VV=F(NN,XX);
% DVV=DF(NN,XX);
% DDVV=DDF(NN,XX);
% 
% DVV0=contour(NN,XX,DVV,[0 0]);

% -- create surf and derivatives
[TT,XX]=meshgrid(t,X);

F = @(t,x) SF(t).*(AF(t).*x + BF(t).*x.^2 + c*(x-CF_(t)).^4 + d*x.^6);
DF = @(t,x) SF(t).*(AF(t) + 2*BF(t).*x + 4*c*(x-CF_(t)).^3 + 6*d*x.^5);
DDF = @(t,x) SF(t).*(2*BF(t) + 12*c*(x-CF_(t)).^2 + 30*d*x.^4);
VV=F(TT,XX);
DVV=DF(TT,XX);
DDVV=DDF(TT,XX);

DVV0=contour(TT,XX,DVV,[0 0]);
DVV0=DVV0(:,2:end);
%%
figure('position',[367   491   938   438]);
surf(NN,XX,VV,'edgecolor','none','facecolor','interp','facealpha',0.9)
colormap parula%jet
zlim(yl)
ylim(xl)
caxis(yl)
hl=camlight('left');
view(-60,45)
xlabel('')
hold on
plot3(DVV0(1,:),DVV0(2,:),F(round(DVV0(1,:)),round(DVV0(2,:))),'k')

%% ===== landscape shift
yl = [-30 15];

figure('position',[367   409   791   520]);

surf(TT,XX,VV,'edgecolor','none','facecolor','interp','facealpha',0.9)
colormap parula%jet
zlim(yl)
ylim(xl)
caxis(yl)
hl=camlight('left');
view(-70,25)
pbaspect([2 1.5 1]);
axis vis3d

hold on
F_0 = F(DVV0(1,:),DVV0(2,:));
% plot3(DVV0(1,:),DVV0(2,:),F_0,'--k')
% -- attractors
attidx=DDF(DVV0(1,:),DVV0(2,:))>0;
DVV0_att = DVV0; DVV0_att(:,~attidx)=nan;
F_0_att = F_0; F_0_att(~attidx)=nan;
plot3(DVV0_att(1,:),DVV0_att(2,:),F_0_att,'--r','linewidth',5)
% -- repellors
DVV0_rpl = DVV0; DVV0_rpl(:,attidx)=nan;
F_0_rpl = F_0; F_0_rpl(attidx)=nan;
plot3(DVV0_rpl(1,:),DVV0_rpl(2,:),F_0_rpl,'--k','linewidth',5)
xlabel('control parameter')
ylabel('state')
zlabel('potential')

% -- plot projection
zplane = -30;
plot3(DVV0_att(1,:),DVV0_att(2,:),repmat(zplane,1,size(DVV0,2)),'r','linewidth',5)
plot3(DVV0_rpl(1,:),DVV0_rpl(2,:),repmat(zplane,1,size(DVV0,2)),'k','linewidth',5)
set(gca,'xtick','','ytick','','ztick','')
box on

% -- plot slice
hSlice=patch(zeros(1,4),[xl fliplr(xl)], repelem(yl,1,2),...
    [0 0 0],'facealpha',0.2,'edgealpha',0,'linewidth',1);

% -- plot ball
r=daspect*0.15;%[range(t) range(X) range(yl)]*0.1;
eps = 0.01;
[xx,yy,zz]=ball(0, 0, F(0,0)+r(3)+eps, r);
hl=surf(xx,yy,zz,'FaceColor','r','EdgeColor','none');

set(gca,'linewidth',5)
%% ===== bifurcation diagram ===== %%
figure
pbaspect([2 1.5 1])
hold on
plot(DVV0_att(1,:),DVV0_att(2,:),'r','linewidth',5)
plot(DVV0_rpl(1,:),DVV0_rpl(2,:),'k','linewidth',5)
set(gca,'xtick','','ytick','','ztick','')
set(gca,'linewidth',2)
xlabel('control parameter')
ylabel('state')
%% ===== landscape ===== %%
figure('position',[367   409   791   520]);
surf(zeros(2,size(X,2))',repmat(X,2,1)',[F(t,X);F(t,X)-1-abs(DF(t,X)*diff(X(1:2)))]','edgecolor','none','facecolor','interp')
colormap parula
hold on
zlim(yl)
ylim(xl)
caxis(yl)
hl=camlight('left');
pbaspect([2 1.5 1]);
axis vis3d




% -- plot ball
r=daspect*0.15;%[range(t) range(X) range(yl)]*0.1;
eps = 0.01;
[xx,yy,zz]=ball(0, 0, F(0,0)+r(3)+eps, r);
hl=surf(xx,yy,zz,'FaceColor','r','EdgeColor','none');

view(-90,0)
xlabel('control parameter')
ylabel('state')
zlabel('potential')
set(gca,'xtick','','ytick','','ztick','')
set(gca,'LineWidth',2)