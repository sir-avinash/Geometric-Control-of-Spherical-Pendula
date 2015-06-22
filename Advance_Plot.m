clear all;
%Lab computer data
dof=[5     10     20   30    40    50    80    100   150   200   250   300   350   400   450   500   550   600   650   700];
time=[0.175 0.211 0.291 0.365 0.444 0.521 0.748 0.904 1.286 1.670 2.054 2.432 2.823 3.198 3.578 3.990 4.350 4.744 5.117 5.511];
min=1/(60*1000);
step=900000*min;
time =time*step;
count=zeros(length(dof)-1,1);
timef=zeros(length(dof)-1,1);
timea=zeros(length(dof)-1,1);
timel=zeros(length(dof)-1,1);
for i=2:20
n1=0;n2=dof(i)/2;n3=0;
n=dof(i);
count(i-1)=142*n-124+139*n-120;
timef(i-1)=(199*n-198+174*n-173+5303.9)/36621;
timea(i-1)=(173*n-229+165*n-256+5303.9)/36621;
mpro(i-1)=(((142*n1+250*n2+358*n3)-124+(139*n1+219*n2+299*n3)-120)+5303.9)/36621;
end
timef =timef*step;
timea =timea*step;
timel =timel*step;
mpro=mpro*step;

set(0,'DefaultLineLineWidth',1.5,'DefaultAxesColorOrder',[0 0 0],...
   'DefaultAxesLineStyleOrder','-|--|:|-.', 'DefaultLineMarkerSize',1.5)

fh1=figure(3);
set(fh1, 'color', 'white'); % sets the color to white 
plot(dof(2:20),mpro,dof(2:20),timef,dof(2:20),timea);
axis([-50 700 -1 130])
set (gca,'fontsize',8,'fontweight','n','fontname','times new romans','linewidth',0.5,'Box', 'off','TickDir','out' );
% title(' Computational Efficiency ','FontSize',8);
xlabel('DOF (n)','FontSize',8);
ylabel('CPU time (min)','FontSize',8);
h=legend('Proposed \theta_1','Fetherstone \theta^2','Mohan and Saha \theta^2_1');
set(h,'Orientation','v','Color', 'none','Box', 'off','Location','best','fontsize',8,'fontweight','n','fontname','times new romans','linewidth',0.5, 'Location', 'SouthEast')
