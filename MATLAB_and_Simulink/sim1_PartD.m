clear
%% parameters
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.81;
%% A,B,C,D
A = zeros(6,6);
A(1,2) = 1;
A(2,3) = m1*g/M;
A(2,5) = m2*g/M;
A(3,4) = 1;
A(4,3) = -(M*g+m1*g)/M/l1;
A(4,5) = -m2*g/M/l1;
A(5,6) = 1;
A(6,3) = -m1*g/M/l2;
A(6,5) = -(M*g+m2*g)/M/l2;
A

B = [0,1/M,0,-1/M/l1,0,-1/M/l2]';
B

C = zeros(3,6); % output x, theta1, theta2
C(1,1) = 1;
C(2,3) = 1;
C(3,5) = 1;
C

D = zeros(3,1);

%% calculate LQR Feedback K
% specify Q and R according to Bryson's rule
xmax = 10000;
xdmax = 100;
theta1_max = 5/180*pi;
theta1_dmax = 100/180*pi;
theta2_max = 5/180*pi;
theta2_dmax = 100/180*pi;
Q = diag([1/xmax^2,1/xdmax^2, ...
    1/theta1_max^2,1/theta1_dmax^2, ...
    1/theta2_max^2,1/theta2_dmax^2]);

u_max = 500;
R = 1/u_max^2;

[K,P,E] = lqr(A,B,Q,R);
K = -K % flip sign
E
%% run linearized and nonlinear SS model simulation
figure(1);clf
set(gcf,'Position',[200,100,1000,1000])
for isim = 1:2
    % set initial condition for simulink model
    x0 = 0;
    xd0 = 0;
    if isim==1 % change initial condition 
        theta1_d0 = 0;
        theta2_d0 = 0;
    else
        theta1_d0 = -1;
        theta2_d0 = 1;
    end
    theta1_0 = 30/180*pi;
    theta2_0 = -30/180*pi;
    
    init_state = [x0,xd0,theta1_0,theta1_d0,theta2_0,theta2_d0]';

    param.M = M;
    param.m1 = m1;
    param.m2 = m2;
    param.l1 = l1;
    param.l2 = l2;
    param.g = g;

    simOut = cell(2,1);
    for imdl = 1:2
        if imdl == 1
            mdl_name = 'linearizedModel';
        else
            mdl_name = 'nonlinearModel';
        end
        simOut{imdl} = sim(mdl_name,'SimulationMode','normal',...
            'SaveState','on','StateSaveName','xout',...
            'SaveOutput','on',...
            'SaveFormat', 'Dataset');
    end
    % plot the difference between linear and nonlinear model output
    color_seq = [1,2,5];
    t1 = simOut{1}.logsout{1}.Values.Time;
    t2 = simOut{2}.logsout{1}.Values.Time;
    for i = 1:3
        subplot(3,2,2*i-1 + (isim-1) )
        co = get(gca,'ColorOrder')
        hold on;
        data1 = simOut{1}.logsout{1}.Values.Data(:,i);
        data2 = simOut{2}.logsout{1}.Values.Data(:,i);
        if i>=2
            data1 = data1/pi*180;
            data2 = data2/pi*180;
        end
        plot(t1,data1,'LineWidth',1,'Color',co(color_seq(i),:),'LineStyle','-');
        plot(t2,data2,'LineWidth',1,'Color',co(color_seq(i),:),'LineStyle','-.');
        legend({'linear','nonlinear'},'Box','off','Location','southeast')
        if isim==2 && i==1
            legend({'linear','nonlinear'},'Box','off','Location','northeast','FontSize',12)
        else
            legend({'linear','nonlinear'},'Box','off','Location','southeast','FontSize',12)
        end
        ylim([1.1*min([data1;data2]),1.1*max([data1;data2])]);
        set(gca,'FontSize',12)
        switch i
            case 1
                title('X','FontSize',16)
                ylabel('Meter')
            case 2
                title('\theta_1','FontSize',16)
                ylabel('Degrees')
                %set(gca,'YTick',abs(theta1_0)/pi*180*[-1,0,1])
            case 3
                title('\theta_2','FontSize',16)
                ylabel('Degrees')
                %set(gca,'YTick',abs(theta2_0)/pi*180*[-1,0,1])
        end
        xlabel('Time (sec)')
        
    end
end
%% save figure
if xmax==10000
    print(gcf,'-dtiff','PartD_fig1');
elseif xmax==100
    print(gcf,'-dtiff','PartD_fig2');
end
%% pring eigenvalues
eig(A+B*K)






