clear
%% set up parameters
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

Bk = [0,1/M,0,-1/M/l1,0,-1/M/l2]';
Bk

% pick output variables
isim = 4
CC = eye(6);
switch isim
    case 1 % output x
        var_sel = 1;
    case 2 % output x, theta1
        var_sel = [1,3];
    case 3 % output x, theta2
        var_sel = [1,5];
    case 4 % output x, theta1, theta2
        var_sel = [1,3,5];
    case 5
        var_sel = [1,2,3,4,5,6];
end
C = CC(var_sel,:);
n_output_var = size(C,1);
C

D = zeros(n_output_var,1);
%% selection matrix S
II = eye(6);
sel_integral = [1];
nc = length(sel_integral);
S = II(sel_integral,:);
%% construct A_bar, B_bar, Q_bar, R_bar
A_bar = [A,zeros(6,nc);S,zeros(nc,nc)];
B_bar = [Bk; zeros(nc,1)];

% Q_bar
xmax = 10000;
xdmax = 100;
theta1_max = 5/180*pi;
theta1_dmax = 100/180*pi;
theta2_max = 5/180*pi;
theta2_dmax = 100/180*pi;
Q = diag([1/xmax^2,1/xdmax^2, ...
    1/theta1_max^2,1/theta1_dmax^2, ...
    1/theta2_max^2,1/theta2_dmax^2]);
Q1=Q;
val = [0.1,10,10,10,10,10];
val = val(sel_integral);
Q2 = diag(val);
Q_bar = blkdiag(Q1,Q2);

% R_bar
u_max = 500;
R = 1/u_max^2;
R_bar = R;

%%  specify noise param
scale_list_p = 1;
scale_list_m = 1;

pow_scale_p = scale_list_p(1)*0.01;
Bd = eye(6);
wn_Ts = 0.01;
wn_pow_x = pow_scale_p*0.001; % process noise for x
wn_pow_theta = pow_scale_p*0.0001; % process noise for theta1 and theta2, 10 fold smaller than x
wn_pow_process_vec = [wn_pow_x;wn_pow_x;...
    wn_pow_theta; wn_pow_theta;...
    wn_pow_theta; wn_pow_theta;];
wn_process_seed = [241;19221;231;313;41;22221];
wn_cov_process = wn_pow_process_vec/wn_Ts;
sigma_d = diag(wn_cov_process);

% specify param for measurement noise
pow_scale_m = scale_list_m(1)*0.1;
wn_pow_measurement_x = pow_scale_m *0.001;
wn_pow_measurement_theta = pow_scale_m*0.0001;
wn_pow_measurement_vec = [wn_pow_measurement_x;wn_pow_measurement_x;
    wn_pow_measurement_theta;wn_pow_measurement_theta;
    wn_pow_measurement_theta;wn_pow_measurement_theta];
wn_pow_measurement_vec = wn_pow_measurement_vec(var_sel);
wn_measurement_seed = [20220;923;55;8887;12132;18];
wn_measurement_seed = wn_measurement_seed(1:n_output_var);
wn_cov_measurement = wn_pow_measurement_vec/wn_Ts;
sigma_v = diag(wn_cov_measurement);
sigma_v_inv = diag(1./wn_cov_measurement);

%% param for the nonlinear model block
param.M = M;
param.m1 = m1;
param.m2 = m2;
param.l1 = l1;
param.l2 = l2;
param.g = g;

%% const reference
xd = [10,0,0,0,0,0]';
u_inf = 0;
%% const disturbance
u_disturb_list = [0,-500];

%% solve for K_bar
[K_bar,P,E] = lqr(A_bar,B_bar,Q_bar,R_bar);
K_bar=-K_bar;
E
K1 = K_bar(1:6);
K2 = K_bar(7:end);
%% solve for L
P = icare(A',C',Bd*sigma_d*Bd',sigma_v,[],[],[]);
L = P*C'*sigma_v_inv;
L
%% simulation
init_state = [0,0,0,0,0,0]';
simOut = cell(1,length(u_disturb_list));
for i = 1:length(u_disturb_list)
    u_disturb = u_disturb_list(i);
    simOut{i} = sim('outputFeedbackCtrl_constRef','SimulationMode','normal',...
        'SaveState','on','StateSaveName','xout',...
        'SaveOutput','on',...
        'SaveFormat', 'Dataset');
end
%% plot result
var_names = {'u','x','\theta_1','\theta_2'};
figure(1)
set(gcf,'Position',[200,100,1600,1200])
clf
plot_data(simOut,var_names,u_disturb_list(2),xd);
subplot(2,2,1)
ylim([-2.5,4.6])
subplot(2,2,2)
ylim([-0.1,18]);
subplot(2,2,3);
ylim([-28,28])
subplot(2,2,4);
ylim([-28,37])

%% save plots
print(gcf, '-dtiff', sprintf("PartG_figConstRef"));

%% calculate tracking error
tracking_err = S*(-inv(A-L*C)*Bk*u_disturb_list(2))


%% helper function
function plot_data(simOut,var_names,u,xd)
    var_sel = [1,3,5];
    color_seq = [4,1,2,5];
    nvar = length(var_sel)+1;
    ncond = length(simOut);
    for ivar = 1:nvar 
        ax(ivar)=subplot(2,2,ivar);
        hold on;
        tend = 70;
        switch ivar
            case 1
                line([0,tend],0*[1,1],'Color','k','LineStyle',':')
                line([0,tend],-u/1000*[1,1],'Color','k','LineStyle',':')
            case 2
                line([0,tend],xd(1)*[1,1],'Color','k','LineStyle',':')
            case 3
                line([0,tend],0*[1,1],'Color','k','LineStyle',':')
            case 4
                line([0,tend],0*[1,1],'Color','k','LineStyle',':')
        end
        for icond = 1:ncond
            % t = simOut{icond}.logsout{2}.Values.Time;
            [x,t1] = find_signal(simOut{icond},'x',var_sel);
            [xhat,t2] = find_signal(simOut{icond},'xhat',var_sel);
            [u,tu] = find_signal(simOut{icond},'u',[]);
            u = u/1000;
        
            if ivar>1
                d1 = x(:,ivar-1);
                d2 = xhat(:,ivar-1);
            else
            end
            if ivar>=3
                % convert to deg
                d1 = d1/pi*180;
                d2 = d2/pi*180;
            end
            
            co = get(gca,'ColorOrder');
            if ivar>1
                if icond==1
                    HP(icond) = plot(t1,d1,'Color',co(color_seq(ivar),:),'LineStyle','-','LineWidth',1.5);
                else
                    HP(icond) = plot(t1,d1,'Color',co(color_seq(ivar),:),'LineStyle','--','LineWidth',1.5);
                end
                % plot(t2,d2,'Color',co(color_seq(ivar),:),'LineStyle','--','LineWidth',1)
            else
                if icond==1
                    HP(icond) = plot(tu,u,'Color',co(color_seq(ivar),:),'LineStyle','-','LineWidth',1.5);
                else
                    HP(icond) = plot(tu,u,'Color',co(color_seq(ivar),:),'LineStyle','--','LineWidth',1.5);
                end

            end
            xlabel('Time (sec)')
            set(gca,'FontSize',16)
           
            switch ivar
                case 1
                    ylabel('Control Effort (\times 1000)')
                case 2
                    ylabel('Meter')
                case 3
                    ylabel('Degrees')
                case 4
                    ylabel('Degrees')
            end
            xlim([0,t1(end)])
        end

        ylim(1.1*ylim);
        legend(HP,{sprintf("$%s$",var_names{ivar}),sprintf("$%s$, with disturbance",var_names{ivar})},...
                    'box','off','Interpreter','latex','FontSize',12,'Location','northeast','FontSize',16)

    end

end

function [data,t] = find_signal(simOut,name,var_sel)
    for i = 1:numElements(simOut.logsout)
        if strcmp(simOut.logsout{i}.Name,name)
            if ~isempty(var_sel)
                data = simOut.logsout{i}.Values.Data;
                data = squeeze(data);
                if size(data,2)~=6
                    data=data';
                end
                data = data(:,var_sel);
            else
                data = simOut.logsout{i}.Values.Data;
                data = data(:);
            end
            t = simOut.logsout{i}.Values.Time;
            return;
        end
    end
end
