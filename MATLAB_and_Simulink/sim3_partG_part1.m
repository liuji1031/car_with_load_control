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

scale_list_p =[1,1,1;
               1,10,100];
scale_list_m = [1,10,100;
                1,1,1];
%%
simOut = cell(size(scale_list_p));
for icond=1:size(scale_list_p,1)
    for iscale = 1:size(scale_list_p,2)
        
        % specify param for uncontrolled input or process noise
        pow_scale_p = scale_list_p(icond,iscale)*0.01;
        Bd = eye(6);
        wn_Ts = 0.01;
        wn_pow_x = pow_scale_p*0.001; % process noise for x
        wn_pow_theta = pow_scale_p*0.0001; % process noise for theta1 and theta2, 10 fold smaller than x
        wn_pow_process_vec = [wn_pow_x;wn_pow_x;...
            wn_pow_theta; wn_pow_theta;...
            wn_pow_theta; wn_pow_theta;];
        wn_process_seed = [23341;191;23111;31313;41;511];
        wn_cov_process = wn_pow_process_vec/wn_Ts;
        sigma_d = diag(wn_cov_process);

        
        % specify param for measurement noise
        pow_scale_m = scale_list_m(icond,iscale)*0.1;
        wn_pow_measurement_x = pow_scale_m *0.001;
        wn_pow_measurement_theta = pow_scale_m*0.0001;
        wn_pow_measurement_vec = [wn_pow_measurement_x;wn_pow_measurement_x;
            wn_pow_measurement_theta;wn_pow_measurement_theta;
            wn_pow_measurement_theta;wn_pow_measurement_theta];
        wn_pow_measurement_vec = wn_pow_measurement_vec(var_sel);
        wn_measurement_seed = [2000;19923;87655;87;1452;1208];
        wn_measurement_seed = wn_measurement_seed(1:n_output_var);
        wn_cov_measurement = wn_pow_measurement_vec/wn_Ts;
        sigma_v = diag(wn_cov_measurement);
        sigma_v_inv = diag(1./wn_cov_measurement);

        % param for the nonlinear model block
        param.M = M;
        param.m1 = m1;
        param.m2 = m2;
        param.l1 = l1;
        param.l2 = l2;
        param.g = g;
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

        [K,P,E] = lqr(A,Bk,Q,R);
        K = -K % flip sign
        E
        %% find L using Kalman-Bucy filter
        % solve the following ricatti equation
        P = icare(A',C',Bd*sigma_d*Bd',sigma_v,[],[],[]);
        L = P*C'*sigma_v_inv;
        L
        eig(A-L*C)
        init_state = [5,-0.5,20/180*pi,-0.5,-20/180*pi,0.5]';
        %init_state = [0,0,20/180*pi,-0.5,-20/180*pi,0.5]';

        %%
        simOut{icond,iscale} = sim('outputFeedbackCtrl','SimulationMode','normal',...
            'SaveState','on','StateSaveName','xout',...
            'SaveOutput','on',...
            'SaveFormat', 'Dataset');
    end
end

%% plot result
var_names = {'u','x','\theta_1','\theta_2'};
for icond=1:size(scale_list_p,1)
    figure(icond)
    set(gcf,'Position',[200,100,1200,1200])
    clf
    plot_data(simOut(icond,:),var_names);
end
%% save plots
for icond=1:size(scale_list_p,1)
    figure(icond)
    print(gcf, '-dtiff', sprintf("PartG_fig%d",icond));
end

%% helper function
function plot_data(simOut,var_names)
    var_sel = [1,3,5];

    color_seq = [4,1,2,5];
    nvar = length(var_sel)+1;
    ncond = length(simOut);
    for icond = 1:ncond
        % t = simOut{icond}.logsout{2}.Values.Time;
        [x,t1] = find_signal(simOut{icond},'x',var_sel);
        [xhat,t2] = find_signal(simOut{icond},'xhat',var_sel);
        [u,tu] = find_signal(simOut{icond},'u',[]);
        for ivar = 1:nvar 
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
            ax(icond,ivar)=subplot(nvar,ncond,icond+(ivar-1)*ncond);
            icond+(ivar-1)*ncond
            hold on;
            co = get(gca,'ColorOrder');
            if ivar>1
                plot(t1,d1,'Color',co(color_seq(ivar),:),'LineStyle','-','LineWidth',1)
                plot(t2,d2,'Color',co(color_seq(ivar),:),'LineStyle','--','LineWidth',1)
            else
                plot(tu,u,'Color',co(color_seq(ivar),:),'LineStyle','-','LineWidth',1)
            end
            xlabel('Time (sec)')
            set(gca,'FontSize',12)

            if ivar>1
                legend({sprintf("$%s$",var_names{ivar}),sprintf('$$\\hat{%s}$$',var_names{ivar})},...
                    'box','off','Interpreter','latex','FontSize',12,'Location','northeast')
            else
                legend({sprintf("$%s$",var_names{ivar})},...
                    'box','off','Interpreter','latex','FontSize',12,'Location','northeast')
            end
            
            switch ivar
                case 1
                    ylabel('Control Effort')
                case 2
                    ylabel('Meter')
                case 3
                    ylabel('Degrees')
                case 4
                    ylabel('Degrees')
            end
        end
    end

    % fix ylim
    for ivar = 1:nvar 
        for icond = 1:ncond
            axes(ax(icond,ivar));
            if icond == 1
                YL = ylim;
            else
                yl_ = ylim;
                YL = [min(YL(1),yl_(1)),max(YL(2),yl_(2))];
            end
        end
        for icond = 1:ncond
            axes(ax(icond,ivar));
            ylim(YL*1.1);
        end
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
