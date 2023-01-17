clear
%% parameters
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.81;

% set up A, B
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

q_hat = 1000;
Q_hat = q_hat*eye(6);
r_hat = 1000;

u_gain = 100;

param.M = M;
param.m1 = m1;
param.m2 = m2;
param.l1 = l1;
param.l2 = l2;
param.g = g;
%% run simulation
simOut = cell(0);
var_names = {'x','\theta_1','\theta_2'};
init_state = [5,-0.5,20/180*pi,-0.5,-20/180*pi,0.5]';
for isim = 1:3
    %% set up parameters
    switch isim
        case 1 % output x
            C = [1,0,0,0,0,0];
            var_sel = 1;
        case 2 % output x, theta2
            C = [1,0,0,0,0,0;
                 0,0,0,0,1,0];
            var_sel = [1,5];
        case 3 % output x, theta1, theta2
            C = [1,0,0,0,0,0;
                 0,0,1,0,0,0;
                 0,0,0,0,1,0];
            var_sel = [1,3,5];
    end
    
    D = zeros(size(C,1),1);
    R_hat = r_hat*eye(size(C,1));
    %% find observer L
    [LT,P,E] = lqr(A',C',Q_hat,R_hat);
    L = LT'
    E
    %% sim with linear model
    simOut{isim,1} = sim('linearizedModel_withObs','SimulationMode','normal',...
        'SaveState','on','StateSaveName','xout',...
        'SaveOutput','on',...
        'SaveFormat', 'Dataset');
    %% sim with nonlinear model
    simOut{isim,2} = sim('nonlinearModel_withObs','SimulationMode','normal',...
        'SaveState','on','StateSaveName','xout',...
        'SaveOutput','on',...
        'SaveFormat', 'Dataset');
    %% plot results
    figure(isim)
    set(gcf,'Position',[200,100,1000,1000])
    clf
    [t_lm,x_lm,xhat_lm,t_nlm,x_nlm,xhat_nlm] = ...
        plot_data(simOut{isim,1},simOut{isim,2},var_names);
    print(gcf, '-dtiff', sprintf("PartF_fig%d",isim));
end
%% helper function
function [t_lm,x_lm,xhat_lm,t_nlm,x_nlm,xhat_nlm] = plot_data(LM,NLM,var_names)
    var_sel = [1,3,5];
    t_lm = LM.logsout{2}.Values.Time;
    x_lm = find_signal(LM,'x',var_sel);
    xhat_lm = find_signal(LM,'xhat',var_sel);

    t_nlm = NLM.logsout{2}.Values.Time;
    x_nlm = find_signal(NLM,'x',var_sel);
    xhat_nlm = find_signal(NLM,'xhat',var_sel);

    color_seq = [1,2,5];
    nvar = length(var_sel);
    for ivar = 1:nvar % two columns, left linear, right nonlinear
        for iplot = 1:2
            if iplot==1
                t = t_lm;
                d1 = x_lm(:,ivar);
                d2 = xhat_lm(:,ivar);
            else
                t = t_nlm;
                d1 = x_nlm(:,ivar);
                d2 = xhat_nlm(:,ivar);
            end
            if ivar>1
                % convert to deg
                d1 = d1/pi*180;
                d2 = d2/pi*180;
            end
            ax(ivar,iplot)=subplot(nvar,2,2*(ivar-1)+iplot);
            
            hold on;
            iname = (var_sel(ivar)+1)/2;
            co = get(gca,'ColorOrder');
            plot(t,d1,'Color',co(color_seq(iname),:),'LineStyle','-','LineWidth',1.5)
            plot(t,d2,'Color',co(color_seq(iname),:),'LineStyle','--','LineWidth',1.5)
            xlabel('Time (sec)')
            set(gca,'FontSize',12)
            if iplot==1
                YL = [min([d1;d2]),max([d1;d2])];
                title(sprintf('Linear system, %s',var_names{ivar}),'FontSize',16)
            else
                YL2 = [min([d1;d2]),max([d1;d2])];
                YL = [min([YL,YL2]),max([YL,YL2])];
                YL = YL * 1.2;
                title(sprintf('Nonlinear system, %s',var_names{ivar}),'FontSize',16)
            end
            
            legend({sprintf("$%s$",var_names{iname}),sprintf('$$\\hat{%s}$$',var_names{iname})},...
                'box','off','Interpreter','latex','FontSize',12,'Location','northeastoutside')
            
            switch ivar
                case 1
                    ylabel('Meter')
                case 2
                    ylabel('Degrees')
                case 3
                    ylabel('Degrees')
            end
        end
        subplot(nvar,2,2*ivar-1);
        ylim(YL);
        subplot(nvar,2,2*ivar);
        ylim(YL);
    end
%     for ivar =  1:nvar
%         pos = get(ax(ivar,2),'Position');
%         pos(1) = pos(1)-0.1;
%         set(ax(ivar,2),'Position',pos);
%     end

end

function data = find_signal(simOut,name,var_sel)
    for i = 1:numElements(simOut.logsout)
        if strcmp(simOut.logsout{i}.Name,name)
            data = simOut.logsout{i}.Values.Data;
            data = squeeze(data);
            if size(data,2)~=6
                data=data';
            end
            data = data(:,var_sel);
            return;
        end
    end
end





