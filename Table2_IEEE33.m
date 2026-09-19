clear; clc; n=32; baseMVA=10;
P_base_MW=[0.1;-0.09;0.12;0.06;0.06;0.2;0.2;0.06;0.06;0.045;-0.06;0.06;0.12;0.06;0.06;0.06;0.09;0.09;0.09;0.09;-0.09;0.09;0.42;0.42;0.06;0.06;0.06;0.12;0.2;0.15;0.21;0.06];
Q_base_Mvar=[0.06;0.04;0.08;0.03;0.02;-0.1;0.1;0.02;0.02;0.03;0.035;0.035;0.08;0.01;0.02;0.02;0.04;0.04;0.04;-0.04;0.04;0.05;0.2;0.2;0.025;0.025;0.02;0.07;0.6;0.07;0.1;0.04];
P_base=P_base_MW/baseMVA; Q_base=Q_base_Mvar/baseMVA; SL0=P_base+1j*Q_base;
%% Ybus
Y = zeros(33,33);
Y(1,2) = -137.9797 + 70.3367j;
Y(2,1) = -137.9797 + 70.3367j;
Y(2,3) = -25.8137 + 13.1477j;
Y(3,2) = -25.8137 + 13.1477j;
Y(2,19) = -51.1502 + 48.8110j;
Y(19,2) = -51.1502 + 48.8110j;
Y(3,4) = -34.7721 + 17.7091j;
Y(4,3) = -34.7721 + 17.7091j;
Y(3,23) = -24.2160 + 16.5465j;
Y(23,3) = -24.2160 + 16.5465j;
Y(4,5) = -33.3937 + 17.0079j;
Y(5,4) = -33.3937 + 17.0079j;
Y(5,6) = -11.2134 + 9.6800j;
Y(6,5) = -11.2134 + 9.6800j;
Y(6,7) = -7.1786 + 23.7293j;
Y(7,6) = -7.1786 + 23.7293j;
Y(6,26) = -62.6890 + 31.9312j;
Y(26,6) = -62.6890 + 31.9312j;
Y(7,8) = -20.3113 + 6.7124j;
Y(8,7) = -20.3113 + 6.7124j;
Y(8,9) = -10.2632 + 7.3736j;
Y(9,8) = -10.2632 + 7.3736j;
Y(9,10) = -10.2183 + 7.2428j;
Y(10,9) = -10.2183 + 7.2428j;
Y(10,11) = -73.4905 + 24.2975j;
Y(11,10) = -73.4905 + 24.2975j;
Y(11,12) = -38.5894 + 12.7601j;
Y(12,11) = -38.5894 + 12.7601j;
Y(12,13) = -6.7435 + 5.3057j;
Y(13,12) = -6.7435 + 5.3057j;
Y(13,14) = -10.8296 + 14.2548j;
Y(14,13) = -10.8296 + 14.2548j;
Y(14,15) = -15.1325 + 13.4682j;
Y(15,14) = -15.1325 + 13.4682j;
Y(15,16) = -14.0065 + 10.2285j;
Y(16,15) = -14.0065 + 10.2285j;
Y(16,17) = -4.4685 + 5.9661j;
Y(17,16) = -4.4685 + 5.9661j;
Y(17,18) = -13.5585 + 10.6319j;
Y(18,17) = -13.5585 + 10.6319j;
Y(19,20) = -5.8806 + 5.2988j;
Y(20,19) = -5.8806 + 5.2988j;
Y(20,21) = -16.5507 + 19.3354j;
Y(21,20) = -16.5507 + 19.3354j;
Y(21,22) = -8.2269 + 10.8775j;
Y(22,21) = -8.2269 + 10.8775j;
Y(23,24) = -10.9933 + 8.6808j;
Y(24,23) = -10.9933 + 8.6808j;
Y(24,25) = -11.0948 + 8.6815j;
Y(25,24) = -11.0948 + 8.6815j;
Y(26,27) = -44.7855 + 22.8025j;
Y(27,26) = -44.7855 + 22.8025j;
Y(27,28) = -8.5152 + 7.5077j;
Y(28,27) = -8.5152 + 7.5077j;
Y(28,29) = -11.3305 + 9.8709j;
Y(29,28) = -11.3305 + 9.8709j;
Y(29,30) = -25.0756 + 12.7725j;
Y(30,29) = -25.0756 + 12.7725j;
Y(30,31) = -8.3211 + 8.2238j;
Y(31,30) = -8.3211 + 8.2238j;
Y(31,32) = -21.8863 + 25.5094j;
Y(32,31) = -21.8863 + 25.5094j;
Y(32,33) = -13.7531 + 21.3839j;
Y(33,32) = -13.7531 + 21.3839j;

for i = 1:33
    Y(i,i) = -sum(Y(i,:));
end
Bsh = [0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000; 0.0000];
Y = Y + diag(1j*Bsh);
%% Partition
PQ=[2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
S=1;
YLL=Y(PQ,PQ);
YLS=Y(PQ,S);
%% ====================== Voltage of generator buses and slack bus ======================
V_M = 1;
V_theta = 0;                         
VS = V_M .* exp(1j*V_theta);
%% ====================== Open-circuit voltage ======================
E = -inv(YLL)*YLS*VS;

%% ====================== Criterion 1： ======================
w1 = NaN;
for w = 1:1e-4:10
    SL = w*SL0;
    A = diag(1./E)*inv(YLL)*diag(conj(SL))*diag(1./conj(E));

    if 4*norm(A,'inf') > 1

        w1 = w - 1e-4;
        break;
    end
end



%% ====================== Criterion 2： ======================
w2 = NaN;
for w = 1:1e-4:10
    SL = w*SL0;
    A = diag(1./E)*inv(YLL)*diag(conj(SL))*diag(1./conj(E));

    a = norm(A,'inf');
    r = norm(sum(A,2),'inf');
    p2 = sqrt(a) + sqrt(r);
    if p2 > 1
        w2 = w - 1e-4;
        break;
    end
end
%% ====================== Criterion 3： ======================
w3 = NaN;
for w = 1:1e-4:10
     SL = w*SL0;
    A  = diag(1./E) * inv(YLL) * diag(conj(SL)) * diag(1./conj(E));
    eta = A * ones(n,1);         
    xi  = abs(A) * ones(n,1);     
    eta_m = norm(eta,'inf');      
    xi_m  = norm(A,'inf');        
    r = 2*(xi + real(eta)) - xi.^2 - abs(eta).^2;
    r_m = max(r);
    p4_1 = r_m + 2*xi_m*eta_m;
    p4_2 = xi_m - eta_m;
    if p4_1 >= 1 || p4_2 > 1
        w3 = w - 1e-4;
        break;
    end
end

%% ====================== Criterion 4： ======================
w4 = NaN;
for w = 1:1e-4:5
    SL = w*SL0;
    A  = diag(1./E) * inv(YLL) * diag(conj(SL)) * diag(1./conj(E));

eta = A*ones(n,1);
xi  = abs(A)*ones(n,1);

gamma = 2*(xi+real(eta))-xi.^2-abs(eta).^2;
gamma_m=max(gamma);

xi_m=norm(xi,inf);

if gamma_m>1
    fail
end

r_radius=sqrt((1-gamma_m)/(2*xi_m^2));

cond1 = all(r_radius*xi < abs(ones(n,1)-eta));

Delta = abs(ones(n,1)-eta).^2-r_radius^2*xi.^2;

h = Delta-(ones(n,1)-eta);

cond2_vec = abs(A*(h./Delta)) + abs(A)*(r_radius*xi./Delta);
cond2 = all(cond2_vec < r_radius*xi);


if ~(cond1 && cond2)
    w4=w-1e-4;
    break;
end
end

fprintf('\n========== 4个判据给出的负载系数(33) ==========' );
fprintf('\nw1* = %.4f', w1);
fprintf('\nw2* = %.4f', w2);
fprintf('\nw3* = %.4f', w3);
fprintf('\nw4* = %.4f', w4);

%% ====================== Theorem 4 & Theorem 5 ======================
SL_base=SL0;

Bmat=diag(1./E)*(YLL\eye(n))*diag(1./conj(E));

K_max=20;
w_range=1:1e-4:10;

w_T4=zeros(K_max+1,1);
w_T5=zeros(K_max+1,1);

for type=1:2

    w_record=zeros(K_max+1,1);

    for iw=1:length(w_range)

        w=w_range(iw);

        SL=w*SL_base;

        A=Bmat*diag(conj(SL));
        AbsA=abs(A);

        one=ones(n,1);

        eta=A*one;
        xi=AbsA*one;


        if type==1

            gamma=2*(xi+real(eta))-xi.^2-abs(eta).^2;
            gamma_m=max(gamma);

            if gamma_m>=1
                continue
            end

            r2=sqrt((1-gamma_m)/(2*norm(xi,inf)^2));

            c=one-eta;
            rho=r2*xi;

        else

            r3=(norm(eta,inf)/2)^(1/3);

            c=one;
            rho=r3*one;

        end


        for k=0:K_max

            den=abs(c).^2-rho.^2;

            if any(den<=0)
                break
            end


            Q1=abs(one-c-A*(c./den));
            Q2=AbsA*(rho./den);

            flag=all(abs(c)>rho) && all(Q1+Q2<rho);


            if flag
                w_record(k+1)=max(w_record(k+1),w);
            end


            if k<K_max

                c=one-A*(c./den);

                rho=AbsA*(rho./den);

            end

        end

    end


    if type==1
        w_T4=w_record;
    else
        w_T5=w_record;
    end

end


fprintf('\n========== Theorem 4 ==========\n')
for K=0:K_max
    fprintf('K=%d  w=%.6f\n',K,w_T4(K+1));
end


fprintf('\n========== Theorem 5 ==========\n')
for K=0:K_max
    fprintf('K=%d  w=%.6f\n',K,w_T5(K+1));
end


%% ====================== CPF  ======================

w_cert = max([w_T4(end), w_T5(end)]);

lambda_cpf = NaN;
w_fail = NaN;

x = ones(n,1);

dw = 1e-6;
sigma_tol = 1e-3;

last_sigma_min = NaN;

for w = w_cert:dw:10
    SL = w*SL0;
    A = diag(1./E)*inv(YLL)*diag(conj(SL))*diag(1./conj(E));
    x_iter=x;
    flag=0;

    for k=1:100

        R=x_iter+A*(1./conj(x_iter))-ones(n,1);
        R2=conj(x_iter)+conj(A)*(1./x_iter)-ones(n,1);

        J_aug=[eye(n),-A*diag(1./conj(x_iter).^2);
              -conj(A)*diag(1./x_iter.^2),eye(n)];

        dx=-J_aug\[R;R2];

        if any(~isfinite(dx))
            break
        end

        x_new=x_iter+dx(1:n);

        if norm(R,'inf')<1e-10 && max(abs(dx(1:n)))<1e-10

            x_iter=x_new;
            flag=1;
            break

        end

        x_iter=x_new;

    end

    R_final=x_iter+A*(1./conj(x_iter))-ones(n,1);

    res=norm(R_final,'inf');

    if flag==0 || res>1e-6
        w_fail=w;
        break;
    end

    x=x_iter;

    J_aug=[eye(n),-A*diag(1./conj(x).^2);
          -conj(A)*diag(1./x.^2),eye(n)];

    sv=svd(J_aug);

    sigma_max_current=sv(1);
    sigma_min_current=sv(end);

    last_sigma_min=sigma_min_current;

    if sigma_min_current < sigma_tol

        lambda_cpf=w;
        x_cpf=x;
        SL_cpf=SL;
        A_cpf=A;
        J_cpf=J_aug;

        sigma_min_cpf=sigma_min_current;
        sigma_max_cpf=sigma_max_current;

        break

    end

end

if isnan(lambda_cpf)

    fprintf('\nCPF did not reach singularity criterion.\n');

else

    VL_cpf=E.*x_cpf;

    R_norm=x_cpf+A_cpf*(1./conj(x_cpf))-ones(n,1);

    err_norm=norm(R_norm,'inf');

    R_pf=YLL*VL_cpf+YLS*VS+conj(SL_cpf)./conj(VL_cpf);

    err_pf=norm(R_pf,'inf');

    fprintf('\n========== CPF boundary ==========');
    fprintf('\nlambda_cpf = %.10f',lambda_cpf);
    fprintf('\nsigma_min(J_aug)=%.4e',sigma_min_cpf);
    fprintf('\nsigma_max(J_aug)=%.4e',sigma_max_cpf);
    fprintf('\nnormalized residual=%.4e',err_norm);
    fprintf('\noriginal PFE residual=%.4e\n',err_pf);

end

%% ========== Generate the required variables for the spectral radius plots of the six criteria ==========
% w5 and w6 denote the certified loading factors obtained from Theorems 4 and 5,
% respectively, with the maximum recursion depth K=20.
w5 = w_T4(end);
w6 = w_T5(end);

% CPFboundary
wCPF = lambda_cpf;

if ~isfinite(wCPF) || wCPF <= 0
    error('CPF 未得到有效边界 lambda_cpf，无法绘制谱半径曲线。');
end

% rho_sp(|A|*diag(|x*|^{-2}))
N_plot = 500;
w_plot = linspace(0, wCPF, N_plot);
rho_plot = zeros(size(w_plot));

x_plot = ones(n,1);

for iw_plot = 1:numel(w_plot)
    w_now = w_plot(iw_plot);
    SL_now = w_now * SL0;
    A_now = Bmat * diag(conj(SL_now));
    
    if iw_plot == 1
        x_plot = ones(n,1);
    else
        converged_plot = false;

        for it_plot = 1:100
            R1_plot = x_plot + A_now*(1./conj(x_plot)) - ones(n,1);
            R2_plot = conj(x_plot) + conj(A_now)*(1./x_plot) - ones(n,1);

            J_plot = [eye(n), -A_now*diag(1./conj(x_plot).^2); ...
                     -conj(A_now)*diag(1./x_plot.^2), eye(n)];

            dx_plot = -J_plot\[R1_plot; R2_plot];

            if any(~isfinite(dx_plot))
                break;
            end

            x_new_plot = x_plot + dx_plot(1:n);

            if norm(R1_plot,'inf') < 1e-11 && ...
                    max(abs(dx_plot(1:n))) < 1e-11
                x_plot = x_new_plot;
                converged_plot = true;
                break;
            end

            x_plot = x_new_plot;
        end

        % Check the final residual
        Rcheck_plot = x_plot + A_now*(1./conj(x_plot)) - ones(n,1);
        if norm(Rcheck_plot,'inf') < 1e-8
            converged_plot = true;
        end

        if ~converged_plot
            error('谱半径曲线在 w = %.6f 处潮流迭代未收敛。', w_now);
        end
    end

    M_rho = abs(A_now) * diag(1./abs(x_plot).^2);
    rho_plot(iw_plot) = max(abs(eig(M_rho)));
end

fprintf('\n========== Spectral-radius plotting data ==========');
fprintf('\nw5 (Theorem 4, K=%d) = %.6f', K_max, w5);
fprintf('\nw6 (Theorem 5, K=%d) = %.6f', K_max, w6);
fprintf('\nwCPF = %.6f', wCPF);
fprintf('\nrho_sp at wCPF = %.6f\n', rho_plot(end));

%% ========== Spectral radius plots of the six criteria==========

required_vars_ieee = {'w_plot','rho_plot','lambda_cpf','w1','w2','w3','w4','w5','w6','wCPF'};
for k_req_ieee = 1:numel(required_vars_ieee)
    if exist(required_vars_ieee{k_req_ieee}, 'var') ~= 1
        error('缺少变量 %s。请先运行原计算程序，再执行本重绘脚本。', ...
              required_vars_ieee{k_req_ieee});
    end
end

if ~isvector(w_plot) || ~isvector(rho_plot) || ...
        numel(w_plot) ~= numel(rho_plot) || isempty(w_plot) || ...
        ~isreal(w_plot) || ~isreal(rho_plot) || ...
        any(~isfinite(w_plot(:))) || any(~isfinite(rho_plot(:)))
    error('w_plot 与 rho_plot 必须是长度相同、非空且有限的实数向量。');
end
if ~isscalar(lambda_cpf) || ~isreal(lambda_cpf) || ...
        ~isfinite(lambda_cpf) || lambda_cpf <= 0
    error('lambda_cpf 必须是有限的正实数。');
end

%% ---------------------- Color schemes and line styles ----------------------
color_curve_ieee = [0, 0, 0];               % 黑色

color_w_ieee = [ ...
    1.00, 0.55, 0.75;   % w1  粉色
    1.00, 0.60, 0.25;   % w2  橙色
    0.00, 0.48, 0.05;   % w3  深绿
    0.00, 0.85, 0.90;   % w4  青色
    0.00, 0.95, 0.05;   % w5  亮绿
    0.00, 0.10, 1.00];  % w6  亮蓝

color_cpf_ieee = [0.78, 0.00, 0.78]; % wCPF 紫红

style_w_ieee = {'--', '--', '--', '--', '--', '--'};
width_w_ieee = [1.40, 1.40, 1.40, 1.40, 1.40, 1.40];

fig_ieee = figure('Color','w', ...
                  'Name','Spectral radius - six criteria', ...
                  'NumberTitle','off', ...
                  'Units','centimeters', ...
                  'Position',[3, 3, 8.8, 6.6], ...
                  'Renderer','painters');
ax_ieee = axes('Parent',fig_ieee, ...
               'Units','normalized', ...
               'Position',[0.19, 0.18, 0.77, 0.77]);
hold(ax_ieee, 'on');

h_curve_ieee = plot(ax_ieee, w_plot, rho_plot, ...
                   'Color',color_curve_ieee, ...
                   'LineStyle','-', 'LineWidth',1.55, ...
                   'DisplayName','$\rho_{\mathrm{sp}}$');

set(ax_ieee, 'Color','w', ...
             'FontName','Times New Roman', ...
             'FontSize',9, ...
             'FontWeight','normal', ...
             'XColor',[0.15,0.15,0.15], ...
             'YColor',[0.15,0.15,0.15], ...
             'LineWidth',0.75, ...
             'TickDir','in', ...
             'TickLength',[0.015,0.015], ...
             'Box','on', ...
             'XGrid','off', 'YGrid','on', ...
             'XMinorGrid','off', 'YMinorGrid','off', ...
             'GridColor',[0.86,0.87,0.89], ...
             'GridAlpha',1, ...
             'GridLineStyle','-');
grid on;

xlim(ax_ieee, [2.8,4]);
ylim(ax_ieee, [0.4,1.2]);

ylim_ieee = [0,1.2];

% 谱半径临界值 rho = 1
h_rho1_ieee = yline(ax_ieee, 1, '-', ...
                     'Color',[1.00,0.00,0.00], ...
                     'LineWidth',1.05, ...
                     'HandleVisibility','off');

xlabel(ax_ieee, '$w$', 'Interpreter','latex', 'FontSize',10);
ylabel(ax_ieee, ...
       '$\rho_{\mathrm{sp}}\left(|A|\,\mathrm{diag}(|x^{*}|^{-2})\right)$', ...
       'Interpreter','latex', 'FontSize',10);

%% ------------------Six criterion vertical lines, labels, and vertical annotations ------------------
w_marks_ieee = [w1, w2, w3, w4, w5, w6];
x_pad_ieee = 0.010*diff(xlim(ax_ieee));

xlim_now_ieee = xlim(ax_ieee);
x_rho_label_ieee = max(xlim_now_ieee(1) + 0.05*diff(xlim_now_ieee), w1 - 0.35);
y_rho_label_ieee = 1.03;
text(ax_ieee, x_rho_label_ieee, y_rho_label_ieee, '$\rho=1$', ...
     'Interpreter','latex', ...
     'FontName','Times New Roman', ...
     'FontSize',8.5, ...
     'Color',[1.00,0.00,0.00], ...
     'HorizontalAlignment','left', ...
     'VerticalAlignment','bottom', ...
     'BackgroundColor','w', ...
     'EdgeColor','none', ...
     'Margin',0.5, ...
     'Clipping','on');

label_height_ieee = [0.93, 0.84, 0.75, 0.66, 0.57, 0.48];

h_w_ieee = gobjects(6,1);

for k_mark_ieee = 1:6
    if ~isfinite(w_marks_ieee(k_mark_ieee)) || ...
            ~isreal(w_marks_ieee(k_mark_ieee))
        warning('判据 w%d 未得到有限实数，跳过对应竖线。', k_mark_ieee);
        continue;
    end
    h_w_ieee(k_mark_ieee) = xline(ax_ieee, w_marks_ieee(k_mark_ieee), ...
                        style_w_ieee{k_mark_ieee}, ...
                        'Color',color_w_ieee(k_mark_ieee,:), ...
                        'LineWidth',width_w_ieee(k_mark_ieee), ...
                        'DisplayName',sprintf('$w_%d$',k_mark_ieee), ...
                        'HandleVisibility','on');
    if isprop(h_w_ieee(k_mark_ieee), 'Alpha')
        h_w_ieee(k_mark_ieee).Alpha = 1;
    end
end

h_vline_cpf_ieee = xline(ax_ieee, wCPF, '--', ...
      'Color',color_cpf_ieee, ...
      'LineWidth',1.40, ...
      'DisplayName','$w_{\mathrm{CPF}}$', ...
      'HandleVisibility','on');

uistack(h_curve_ieee, 'top');

legend_handles_ieee = [h_w_ieee; h_vline_cpf_ieee];
legend_labels_ieee = {sprintf('$w_1 = %.4f$', w1), ...
                      sprintf('$w_2 = %.4f$', w2), ...
                      sprintf('$w_3 = %.4f$', w3), ...
                      sprintf('$w_4 = %.4f$', w4), ...
                      sprintf('$w_5 = %.4f$', w5), ...
                      sprintf('$w_6 = %.4f$', w6), ...
                      sprintf('$w_{\\mathrm{CPF}} = %.4f$', wCPF)};

lgd_ieee = legend(ax_ieee, legend_handles_ieee, legend_labels_ieee, ...
       'Interpreter','latex', ...
       'Location','northwest', ...
       'FontName','Times New Roman', ...
       'FontSize',8.0, ...
       'Box','on');
lgd_ieee.LineWidth = 0.6;

hold(ax_ieee, 'off');
drawnow;