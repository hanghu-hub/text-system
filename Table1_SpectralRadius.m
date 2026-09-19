%% Three-bus power-flow random test: independently generate 100000 cases per run with full batch processing
% Requires MATLAB R2023a or later: pageeig and pagemldivide are used.
% No additional toolbox is required. No sample-wise loop or arrayfun/cellfun wrapper is used.
% The only while loop advances Newton iterations and updates all candidate solutions simultaneously.
%
% Buses 1 and 2 are PQ buses; bus 3 is the slack bus with V3=1. No shunt elements are included.
% The three columns of Z correspond to z12, z13, and z23. S=P+1i*Q uses the positive consumed power convention.

% Procedure: batch generation of sixth-order polynomials -> all roots by pageeig -> voltage recovery -> batch refinement.
%      -> residual filtering, duplicate removal -> batch rho_sp calculation and statistics.
% rho_sp = rho(abs(A)*diag(1./abs(V).^2)), A=Y\diag(conj(S))。
% Retain the line-circle intersection recovery branch for D(t)=0 to avoid missing solutions in symmetric cases.
%
% The workspace variable result stores all inputs, solutions, and rho_sp values; a MAT file is saved separately.
% result.V1(k,j)、V2(k,j)、rho_sp(k,j)：第 k 套的第 j 个解。
% Each case is sorted by rho_sp; empty entries are filled with NaN, usually stored as an N x 6 matrix.
% result.status==0: no numerical enumeration issues detected; nonzero values require further verification.
% This is a double-precision numerical all-root test, not a symbolic proof. Unresolved cases are not counted as confirmed infeasible cases.

N = 100000;
R_range = [0.1, 1.0];                 % Resistive-inductive branches: R>0, X>0
X_range = [0.1, 1.0];
P_range = [-0.6, 0.6];                % The two PQ nodes are sampled independently.
Q_range = [-0.6, 0.6];                % For pure load cases, both P and Q ranges can be set to positive values.
seed = [];                            % []: generate new data each run; integer: use a fixed random seed.
save_results = true;

assert(~verLessThan('matlab','9.14'), ...
    '此批处理版本需要 MATLAB R2023a 或更新版本（pageeig）。');
options = struct('Rrange',R_range,'Xrange',X_range, ...
    'Prange',P_range,'Qrange',Q_range,'seed',seed, ...
    'rhoTol',1e-8,'realTol',1e-9,'nearRealTol',1e-5, ...
    'clusterTol',1e-7,'denTol',1e-9,'resTol',1e-10, ...
    'polishTol',5e-14,'newtonMax',8,'duplicateTol',1e-8);

result = pf3_run_random(N,options);

if save_results
    output_folder = fileparts(mfilename('fullpath'));
    if isempty(output_folder)
        output_folder = pwd;
    end
    result.output_file = fullfile(output_folder, ...
        ['pf3_random_results_',datestr(now,'yyyymmdd_HHMMSS_FFF'),'.mat']);
    save_timer = tic;
    save(result.output_file,'result','-v7');
    fprintf('\n结果已保存：%s\n保存耗时：%.3f 秒\n', ...
        result.output_file,toc(save_timer));
end


function result = pf3_run_random(N,opt)
    assert(isscalar(N) && N>=1 && N==fix(N),'N 必须为正整数。');
    assert(numel(opt.Rrange)==2 && all(isfinite(opt.Rrange)) ...
        && opt.Rrange(1)>0 && diff(opt.Rrange)>0);
    assert(numel(opt.Xrange)==2 && all(isfinite(opt.Xrange)) ...
        && opt.Xrange(1)>0 && diff(opt.Xrange)>0);
    assert(numel(opt.Prange)==2 && all(isfinite(opt.Prange)) && diff(opt.Prange)>0);
    assert(numel(opt.Qrange)==2 && all(isfinite(opt.Qrange)) && diff(opt.Qrange)>0);
    pf3_self_check(opt);
    if isempty(opt.seed)
        rng('shuffle');
    else
        rng(opt.seed,'twister');
    end
    rng_state = rng;
    Z = opt.Rrange(1)+diff(opt.Rrange)*rand(N,3) ...
        +1i*(opt.Xrange(1)+diff(opt.Xrange)*rand(N,3));
    S = opt.Prange(1)+diff(opt.Prange)*rand(N,2) ...
        +1i*(opt.Qrange(1)+diff(opt.Qrange)*rand(N,2));

    fprintf('\n开始批量求解 %d 套随机数据。\n',N);
    fprintf('R ∈ [%.3g, %.3g], X ∈ [%.3g, %.3g]\n',opt.Rrange,opt.Xrange);
    fprintf('P ∈ [%.3g, %.3g], Q ∈ [%.3g, %.3g]\n',opt.Prange,opt.Qrange);
    calculation_timer = tic;
    r = pf3_solve_batch(Z,S,opt);

    % Abnormal cases are recomputed in batch: exchange PQ nodes and eliminate using |V2|^2.
    retried = r.status~=0;
    retry_ids = find(retried);
    if ~isempty(retry_ids)
        fprintf('对 %d 套异常样本交换 PQ 节点，批量复核。\n',numel(retry_ids));
        w = pf3_solve_batch(Z(retry_ids,[1,3,2]),S(retry_ids,[2,1]),opt);
        temporary = w.V1;
        w.V1 = w.V2;
        w.V2 = temporary;
        covered = pf3_subset(r.V1(retry_ids,:),r.V2(retry_ids,:), ...
            r.good(retry_ids,:),w.V1,w.V2,w.good,opt.duplicateTol);
        accepted = w.status==0 & covered;
        ids = retry_ids(accepted);
        r.V1(ids,:) = w.V1(accepted,:);
        r.V2(ids,:) = w.V2(accepted,:);
        r.rho(ids,:) = w.rho(accepted,:);
        r.err(ids,:) = w.err(accepted,:);
        r.good(ids,:) = w.good(accepted,:);
        r.status(ids) = uint16(0);

        % When two enumerations cannot mutually confirm results, merge and retain all valid candidate solutions.
        ids = retry_ids(~accepted);
        if ~isempty(ids)
            merged_v1 = [r.V1(ids,:),w.V1(~accepted,:)];
            merged_v2 = [r.V2(ids,:),w.V2(~accepted,:)];
            merged_rho = [r.rho(ids,:),w.rho(~accepted,:)];
            merged_err = [r.err(ids,:),w.err(~accepted,:)];
            merged_good = pf3_deduplicate(merged_v1,merged_v2, ...
                [r.good(ids,:),w.good(~accepted,:)],opt.duplicateTol);
            extra = size(r.V1,2)+(1:size(w.V1,2));
            r.V1(:,extra)=NaN; r.V2(:,extra)=NaN;
            r.rho(:,extra)=NaN; r.err(:,extra)=NaN; r.good(:,extra)=false;
            r.V1(ids,:)=merged_v1; r.V2(ids,:)=merged_v2;
            r.rho(ids,:)=merged_rho; r.err(ids,:)=merged_err;
            r.good(ids,:)=merged_good;
            r.status(ids)=bitor(bitor(r.status(ids),w.status(~accepted)),uint16(16));
        end
    end

    % Sort each row and compress empty entries; all samples are processed simultaneously.
    r = pf3_pack(r);
    elapsed = toc(calculation_timer);
    n_solutions = sum(r.good,2);
    n_rho_le1 = sum(r.good & r.rho<=1,2);
    n_rho_strict = sum(r.good & r.rho<1-opt.rhoTol,2);
    n_rho_boundary = sum(r.good & abs(r.rho-1)<=opt.rhoTol,2);
    n_rho_with_tol = sum(r.good & r.rho<=1+opt.rhoTol,2);
    checked = r.status==0;
    histogram = accumarray(n_solutions(checked)+1,ones(nnz(checked),1),[7,1]);

    result.N = N;
    result.options = opt;
    result.rng_state = rng_state;
    result.z12=Z(:,1); result.z13=Z(:,2); result.z23=Z(:,3);
    result.S1=S(:,1); result.S2=S(:,2);
    result.V1=r.V1; result.V2=r.V2;
    result.rho_sp=r.rho; result.power_residual=r.err;
    result.n_solutions=uint8(n_solutions); result.status=r.status;
    result.used_second_elimination=retried;
    result.n_rho_le1=n_rho_le1;
    result.n_rho_le1_with_tol=n_rho_with_tol;
    result.n_rho_strict=n_rho_strict;
    result.n_rho_boundary=n_rho_boundary;
    result.suspicious_indices=find(n_rho_with_tol>1);
    result.strict_violation_indices=find(n_rho_strict>1);
    result.boundary_indices=find(n_rho_boundary>0);
    result.unresolved_indices=find(~checked);
    result.summary=struct('total',N,'checked',nnz(checked), ...
        'unresolved',nnz(~checked),'solvable',nnz(checked & n_solutions>0), ...
        'solution_count_histogram_0_to_6',histogram, ...
        'total_validated_solutions',sum(n_solutions), ...
        'max_count_rho_le1',max(n_rho_le1), ...
        'max_count_rho_le1_with_tol',max(n_rho_with_tol), ...
        'possible_violations',nnz(n_rho_with_tol>1), ...
        'strict_violations',nnz(n_rho_strict>1),'compute_seconds',elapsed);

    fprintf('\n总套数：%d；数值枚举正常：%d；待复核：%d\n',N,nnz(checked),nnz(~checked));
    fprintf('以下解数分布只包含数值枚举正常的样本：\n');
    disp(table((0:6).',histogram,'VariableNames',{'Solutions','Cases'}));
    fprintf('每套 rho_sp <= 1 的最大解数：%d\n',max(n_rho_le1));
    fprintf('计入 %.1e 容差后的最大解数：%d\n',opt.rhoTol,max(n_rho_with_tol));
    fprintf('出现两个以上合格候选解的套数：%d\n',nnz(n_rho_with_tol>1));
    fprintf('rho_sp 接近 1 的样本数：%d\n',nnz(n_rho_boundary>0));
    fprintf('计算耗时：%.3f 秒（不含 MAT 文件保存）。\n',elapsed);
    if isempty(result.suspicious_indices)
        fprintf('本次已恢复的解中，未发现两个不同解同时满足 rho_sp <= 1。\n');
    else
        fprintf('疑似反例编号保存在 result.suspicious_indices：\n');
        disp(result.suspicious_indices(1:min(end,20)).');
    end
    if any(~checked)
        fprintf('result.unresolved_indices 中的样本仍需复核，不能记为已确认无解。\n');
    end
end

function r = pf3_solve_batch(Z,S,opt)
    N=size(Z,1);
    p=pf3_polynomials(Z,S);
    [roots_t,bad_polynomial]=pf3_page_roots(p.F);
    t=real(roots_t);
    root_scale=max(abs(t),1e-10);
    positive=t>0 & abs(imag(roots_t))<=opt.nearRealTol*root_scale;
    near_real=positive & abs(imag(roots_t))>opt.realTol*root_scale;

    [sorted_t,order]=sort(t,2);
    order_index=(1:N).'+N*(order-1);
    sorted_positive=positive(order_index);
    gaps=abs(diff(sorted_t,1,2))<=opt.clusterTol*max( ...
        max(abs(sorted_t(:,1:5)),abs(sorted_t(:,2:6))),1e-10);
    gaps=gaps & sorted_positive(:,1:5) & sorted_positive(:,2:6);
    clustered=false(N,6);
    clustered(order_index)=[gaps,false(N,1)] | [false(N,1),gaps];

    den=pf3_eval2(p.D,t);
    lv=p.L(:,1)+p.L(:,2).*t;
    mv=p.M(:,1)+p.M(:,2).*t;
    den_ok=abs(den)>opt.denTol*max(abs(lv).^2+abs(mv).^2,realmin);
    regular_v1=pf3_eval3(p.X,t)./den;

    % When D(t)=0, L(t)*V1+M(t)*conj(V1)=T(t) defines a line.
    % Combining with the circle |V1|^2=t, each D root provides at most two candidate voltages.
    [dt,plus_v1,minus_v1,plus_ok,minus_ok]=pf3_common_roots(p);
    covered=any(reshape(plus_ok | minus_ok,N,1,2) ...
        & abs(reshape(t,N,6,1)-reshape(dt,N,1,2)) ...
        <=1e-4*max(abs(reshape(dt,N,1,2)),1e-10),3);
    root_issue=any(positive & ~covered & (~den_ok | near_real | clustered),2);
    V1=[regular_v1,plus_v1,minus_v1];
    anchors=[t,dt,dt];
    wanted=[positive & den_ok,plus_ok,minus_ok];
    y=1./Z; b=y(:,1); a=y(:,2); c=y(:,3);
    V2=((a+b).*V1-a+conj(S(:,1))./conj(V1))./b;
    initial=wanted & isfinite(V1) & isfinite(V2) & abs(V1)>0 & abs(V2)>0;

    % Flatten valid candidate solutions; update Newton steps by solving batched linear systems.
    [V1,V2,err,good]=pf3_polish_batch(V1,V2,initial,anchors,y,S,opt);
    failed_recovery=any(wanted & ~good,2);
    good=pf3_deduplicate(V1,V2,good,opt.duplicateTol);
    rho=pf3_rho_batch(V1,V2,Z,S);
    bad_rho=any(good & ~isfinite(rho),2);
    count=sum(good,2);

    % Status flags: 1 coefficient degeneration/zero power; 2 roots near the real axis, repeated roots, or denominator degeneration;
    % 4 back-substitution/refinement/spectral-radius anomaly; 8 odd number of solutions or more than 6 solutions; 16 unconfirmed after two eliminations.
    status=uint16(double(bad_polynomial | any(S==0,2)) ...
        +2*double(root_issue)+4*double(failed_recovery | bad_rho) ...
        +8*double(mod(count,2)~=0 | count>6));
    r=struct('V1',V1,'V2',V2,'rho',rho,'err',err,'good',good,'status',status);
end

%%% Batch sixth-order polynomial coefficients arranged from constant to sixth-order terms
function p = pf3_polynomials(Z,S)
    y=1./Z;
    scale=max(max(abs(y),[],2),max(abs(S),[],2));
    y=y./scale; S=S./scale;           % 同比例缩放不改变电压解
    b=y(:,1); a=y(:,2); c=y(:,3); s1=S(:,1); s2=S(:,2);
    E=a.*c+a.*b+b.*c;
    b0=s1; b1=conj(a+b);
    h0=(b+c).*conj(s1); h1=E;
    l0=conj(a).*h0; l1=conj(a).*h1;
    m0=E.*b0; m1=E.*b1;
    n0=b0.*h0;
    n1=b0.*h1+b1.*h0+conj(a).*E+abs(b).^2.*conj(s2);
    n2=b1.*h1;
    d0=abs(l0).^2-abs(m0).^2;
    d1=2*real(l0.*conj(l1)-m0.*conj(m1));
    d2=abs(l1).^2-abs(m1).^2;
    x0=n0.*conj(l0)-m0.*conj(n0);
    x1=n0.*conj(l1)+n1.*conj(l0)-m0.*conj(n1)-m1.*conj(n0);
    x2=n1.*conj(l1)+n2.*conj(l0)-m0.*conj(n2)-m1.*conj(n1);
    x3=n2.*conj(l1)-m1.*conj(n2);
    F=[abs(x0).^2, ...
        2*real(x0.*conj(x1))-d0.^2, ...
        abs(x1).^2+2*real(x0.*conj(x2))-2*d0.*d1, ...
        2*real(x0.*conj(x3)+x1.*conj(x2))-d1.^2-2*d0.*d2, ...
        abs(x2).^2+2*real(x1.*conj(x3))-2*d1.*d2, ...
        2*real(x2.*conj(x3))-d2.^2,abs(x3).^2];
    p=struct('F',F,'X',[x0,x1,x2,x3],'D',[d0,d1,d2], ...
        'L',[l0,l1],'M',[m0,m1],'T',[n0,n1,n2]);
end

%%% Compute all roots of N sixth-order polynomials using one pageeig call
function [roots_t,bad] = pf3_page_roots(F)
    N=size(F,1);
    bad=any(~isfinite(F),2) | F(:,7)==0;
    F(bad,:)=repmat([1,0,0,0,0,0,1],nnz(bad),1);
    alpha=min(max(max(abs(F(:,1)./F(:,7)),1e-240).^(1/6),1e-8),1e8);
    ff=F.*alpha.^(0:6);
    top=-ff(:,6:-1:1)./ff(:,7);
    bad=bad | any(~isfinite(top),2);
    top(bad,:)=repmat([0,0,0,0,0,-1],nnz(bad),1);
    C=repmat(diag(ones(5,1),-1),1,1,N);        % 6 x 6 x N
    C(1,:,:)=permute(top,[3,2,1]);
    roots_t=reshape(pageeig(C),6,N).'.*alpha;   % N x 6
    root_error=abs(pf3_eval6(F,roots_t));
    root_size=pf3_eval6(abs(F),abs(roots_t));
    bad=bad | any(~isfinite(roots_t),2) ...
        | any(root_error>1e-8*max(root_size,realmin),2);
end

%%% Batch recovery of solutions when the elimination denominator is zero
function [dt,vp,vm,gp,gm] = pf3_common_roots(p)
    d0=p.D(:,1); d1=p.D(:,2); d2=p.D(:,3);
    discriminant=d1.^2-4*d2.*d0;
    sign_d1=ones(size(d1)); sign_d1(d1<0)=-1;
    q=-0.5*(d1+sign_d1.*sqrt(max(discriminant,0)));
    dt=[q./d2,d0./q];
    common=discriminant>=0 & dt>0 & isfinite(dt);
    common=common & abs(pf3_eval6(p.F,dt)) ...
        <=5e-13*max(pf3_eval6(abs(p.F),abs(dt)),realmin);
    l=p.L(:,1)+p.L(:,2).*dt;
    m=p.M(:,1)+p.M(:,2).*dt;
    n=pf3_eval2(p.T,dt);
    wx1=real(l+m); wy1=-imag(l-m);
    wx2=imag(l+m); wy2=real(l-m);
    choose_first=wx1.^2+wy1.^2>=wx2.^2+wy2.^2;
    wx=wx2; wy=wy2; rhs=imag(n);
    wx(choose_first)=wx1(choose_first);
    wy(choose_first)=wy1(choose_first);
    real_n=real(n); rhs(choose_first)=real_n(choose_first);
    ww=wx.^2+wy.^2;
    cx=rhs.*wx./ww; cy=rhs.*wy./ww;
    height=dt-cx.^2-cy.^2;
    common=common & ww>0 & height>=-1e-10*max(dt,1e-10);
    factor=sqrt(max(height,0)./ww);
    ox=-wy.*factor; oy=wx.*factor;
    vp=cx+ox+1i*(cy+oy); vm=cx-ox+1i*(cy-oy);
    ep=max(abs(wx1.*real(vp)+wy1.*imag(vp)-real(n)), ...
        abs(wx2.*real(vp)+wy2.*imag(vp)-imag(n)));
    em=max(abs(wx1.*real(vm)+wy1.*imag(vm)-real(n)), ...
        abs(wx2.*real(vm)+wy2.*imag(vm)-imag(n)));
    tolerance=1e-8*max(max(abs(n),sqrt(max(dt,0).*ww)),realmin);
    gp=common & ep<=tolerance; gm=common & em<=tolerance;
end

%%% Newton while loop controls only iteration rounds; all active candidates solve 4 x 4 systems simultaneously
function [V1,V2,errors,good] = pf3_polish_batch(V1,V2,initial,t,y,S,opt)
    N=size(V1,1);
    ids=find(initial(:));
    cases=mod(ids-1,N)+1;
    
    u=reshape(V1(ids),[],1); v=reshape(V2(ids),[],1);
    target_t=reshape(t(ids),[],1);
    a=y(cases,2); b=y(cases,1); c=y(cases,3);
    s1=S(cases,1); s2=S(cases,2);
    tolerance_scale=1+max(abs(s1),abs(s2));
    blocked=false(size(ids));
    iteration=0;
    while iteration<opt.newtonMax
        I1=(a+b).*u-b.*v-a; I2=(b+c).*v-b.*u-c;
        g1=conj(u).*I1+conj(s1); g2=conj(v).*I2+conj(s2);
        active=find(max(abs(g1),abs(g2))>opt.polishTol*tolerance_scale & ~blocked);
        if isempty(active)
            break;
        end
        aa=a(active); bb=b(active); cc=c(active);
        uu=u(active); vv=v(active); i1=I1(active); i2=I2(active);
        dx11=i1+conj(uu).*(aa+bb); dx12=-conj(uu).*bb;
        dx21=-conj(vv).*bb; dx22=i2+conj(vv).*(bb+cc);
        dy11=-1i*i1+1i*conj(uu).*(aa+bb); dy12=1i*dx12;
        dy21=1i*dx21; dy22=-1i*i2+1i*conj(vv).*(bb+cc);
        J=reshape([real(dx11),real(dx21),imag(dx11),imag(dx21), ...
            real(dx12),real(dx22),imag(dx12),imag(dx22), ...
            real(dy11),real(dy21),imag(dy11),imag(dy21), ...
            real(dy12),real(dy22),imag(dy12),imag(dy22)].',4,4,[]);
        rhs=reshape([real(g1(active)),real(g2(active)), ...
            imag(g1(active)),imag(g2(active))].',4,1,[]);
        [delta,reciprocal_condition]=pagemldivide(J,rhs);
        delta=reshape(delta,4,[]).';
        safe=reshape(reciprocal_condition,[],1)>1e-14 & all(isfinite(delta),2);
        blocked(active(~safe))=true;
        update=active(safe);
        u(update)=u(update)-delta(safe,1)-1i*delta(safe,3);
        v(update)=v(update)-delta(safe,2)-1i*delta(safe,4);
        iteration=iteration+1;
    end
    V1(ids)=u; V2(ids)=v;
    err=max(abs(conj(u).*((a+b).*u-b.*v-a)+conj(s1)), ...
        abs(conj(v).*((b+c).*v-b.*u-c)+conj(s2)));
    valid=isfinite(u) & isfinite(v) & abs(u)>0 & abs(v)>0 ...
        & err<=opt.resTol*tolerance_scale ...
        & abs(abs(u).^2-target_t)<=1e-5*max(abs(target_t),1e-10);
    good=false(size(V1)); good(ids)=valid;
    errors=nan(size(V1)); errors(ids)=err;
end

%%% Batch duplicate removal: pairwise distance arrays are built only for potentially duplicated rows
function good = pf3_deduplicate(V1,V2,good,tol)
    t=abs(V1).^2; t(~good)=Inf;
    ts=sort(t,2);
    rows=find(any(isfinite(ts(:,2:end)) ...
        & diff(ts,1,2)<1e-7*(1+ts(:,2:end)),2));
    if isempty(rows)
        return;
    end
    M=numel(rows); K=size(V1,2);
    u=V1(rows,:); v=V2(rows,:); ok=good(rows,:);
    u1=reshape(u,M,K,1); u2=reshape(u,M,1,K);
    v1=reshape(v,M,K,1); v2=reshape(v,M,1,K);
    distance=max(abs(u1-u2),abs(v1-v2));
    scale=1+max(max(abs(u1),abs(v1)),max(abs(u2),abs(v2)));
    duplicate=distance<=tol*scale ...
        & reshape(ok,M,K,1) & reshape(ok,M,1,K) ...
        & reshape(tril(true(K),-1),1,K,K);
    good(rows,:)=ok & ~any(duplicate,3);
end

%%% rho_sp: Perron root of the nonnegative 2 x 2 matrix is computed element-wise
function rho = pf3_rho_batch(V1,V2,Z,S)
    y=1./Z; b=y(:,1); a=y(:,2); c=y(:,3);
    E=a.*c+a.*b+b.*c;
    A11=abs((b+c)./E).*abs(S(:,1)); A12=abs(b./E).*abs(S(:,2));
    A21=abs(b./E).*abs(S(:,1)); A22=abs((a+b)./E).*abs(S(:,2));
    m11=A11./abs(V1).^2; m12=A12./abs(V2).^2;
    m21=A21./abs(V1).^2; m22=A22./abs(V2).^2;
    rho=0.5*(m11+m22+sqrt((m11-m22).^2+4*m12.*m21));
end

function covered = pf3_subset(u1,u2,ug,v1,v2,vg,tol)
    M=size(u1,1); K=size(u1,2); L=size(v1,2);
    a=reshape(u1,M,K,1); b=reshape(u2,M,K,1);
    c=reshape(v1,M,1,L); d=reshape(v2,M,1,L);
    distance=max(abs(a-c),abs(b-d));
    scale=1+max(max(abs(a),abs(b)),max(abs(c),abs(d)));
    matched=distance<=tol*scale & reshape(vg,M,1,L);
    covered=all(~ug | any(matched,3),2);
end

function r = pf3_pack(r)
    N=size(r.V1,1);
    key=r.rho;
    key(r.good & ~isfinite(key))=realmax;
    key(~r.good)=Inf;
    [~,order]=sort(key,2);
    width=max(6,max(sum(r.good,2)));
    take=(1:N).'+N*(order(:,1:width)-1);
    r.V1=r.V1(take); r.V2=r.V2(take);
    r.rho=r.rho(take); r.err=r.err(take); r.good=r.good(take);
    r.V1(~r.good)=NaN; r.V2(~r.good)=NaN;
    r.rho(~r.good)=NaN; r.err(~r.good)=NaN;
end

function value = pf3_eval2(c,t)
    value=(c(:,3).*t+c(:,2)).*t+c(:,1);
end

function value = pf3_eval3(c,t)
    value=((c(:,4).*t+c(:,3)).*t+c(:,2)).*t+c(:,1);
end

function value = pf3_eval6(c,t)
    value=c(:,7).*t+c(:,6);
    value=value.*t+c(:,5);
    value=value.*t+c(:,4);
    value=value.*t+c(:,3);
    value=value.*t+c(:,2);
    value=value.*t+c(:,1);
end

function pf3_self_check(opt)
    % Original symmetric six-solution test case, including two non-real voltage solutions from D(t)=0.
    Z=[1+1i,1+1i,1+1i]; S=[-1.5-1.5i,-1.5-1.5i];
    r=pf3_solve_batch(Z,S,opt);
    expected=sort([(7-sqrt(13))/6;(7+sqrt(13))/6;3;3; ...
        (9+sqrt(33))/4;(9+sqrt(33))/4]);
    obtained=sort(r.rho(r.good));
    assert(nnz(r.good)==6,'六解参考案例未能恢复全部六个解。');
    assert(max(abs(obtained(:)-expected))<1e-8,'六解参考案例的 rho_sp 检查失败。');
    % Analytical critical solution V1=V2=1/2 with rho_sp=1; floating-point repeated roots are not treated as different solutions.
    rho_fold=pf3_rho_batch(0.5,0.5,Z,[(1+1i)/8,(1+1i)/8]);
    assert(abs(rho_fold-1)<1e-12,'rho_sp=1 临界参考案例检查失败。');
end
