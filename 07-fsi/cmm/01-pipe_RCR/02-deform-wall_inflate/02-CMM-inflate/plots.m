clear all; close all; clc;

is = 200;
ie = 800;

t = linspace(1,800,800)'.*0.005;

fname = '24-procs/B_CM_Velocity_flux.txt';
fid = fopen(fname,'r');
C = textscan(fid,'%f %f %f','HeaderLines',3);
fclose(fid);
Qi = C{:,1};
Qo = C{:,2};
Qw = C{:,3};

fname = '24-procs/B_CM_Pressure_average.txt';
fid = fopen(fname,'r');
C = textscan(fid,'%f %f %f','HeaderLines',3);
fclose(fid);
Pi = C{:,1};
Po = C{:,2};

t  = t(is:ie);
Qi = Qi(is:ie);
Qo = Qo(is:ie);
Qw = Qw(is:ie);
Pi = Pi(is:ie);
Po = Po(is:ie);

figure;
    set(gcf,'DefaultAxesFontSize',16);
    set(gcf,'DefaultTextFontSize',16);
    plot(t,-Qi,'k-','LineWidth',2); hold on;
    plot(t,-Qw,'b-','LineWidth',2);
    plot(t,Qo,'r--','LineWidth',2); hold off;
    xlabel('t [s]'); ylabel('Flow (ml/s)');
    legend('Inflow', 'Wall flux', 'Outflow');
    
figure;
    set(gcf,'DefaultAxesFontSize',16);
    set(gcf,'DefaultTextFontSize',16);
    plot(t,-(Qi+Qw),'k-','LineWidth',2); hold on;
    plot(t,Qo','r--','LineWidth',2); hold off;
    xlabel('t [s]'); ylabel('Flow (ml/s)');
    legend('Net Inflow', 'Net Outflow');
    
figure;
    set(gcf,'DefaultAxesFontSize',16);
    set(gcf,'DefaultTextFontSize',16);
    plot(t,Pi./1334.0,'k-','LineWidth',2); hold on;
    plot(t,Po./1334.0,'r--','LineWidth',2); hold off;
    xlabel('t [s]'); ylabel('Pressure (mmHg)');
    legend('P in', 'P out');
    
