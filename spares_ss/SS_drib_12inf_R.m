function [i,d_w,choose,t,lamda2] = SS_drib_12inf_R(A, R, s, b0, K0, Q,alpha)
fs = 16000;
K = 256;
Rvv = R{2}; 
M = size(Rvv,1);
%% ATF function construction from the target position
sour = s{1}; itf1 = s{2}; itf2 = s{3};
[att,dly]=steer(sour); 
for k = 1:K/2+1           
    dvec_s(:,k) = att(1:M).*exp(-sqrt(-1)*2*pi*(k-1)/K*fs*dly(1:M)); 
end
% MVDR
zz = 3;
for k = 1:K/2+1           
    Rinv2 = inv( Rvv(zz,zz,k) );
    Wopt2(:,k) = Rinv2 * dvec_s(zz,k) / (dvec_s(zz,k)' * Rinv2 * dvec_s(zz,k) );
    oNoise_max(:,k) = abs( Wopt2(:,k)'*Rvv(zz,zz,k)*Wopt2(:,k) );
end
for k = 1:K/2+1           
    Rinv = inv( Rvv(1:M,1:M,k) );
    Wopt(:,k) = Rinv * dvec_s(:,k) / (dvec_s(:,k)' * Rinv * dvec_s(:,k) );
    oNoise_min(:,k) = abs( Wopt(:,k)'*Rvv(1:M,1:M,k)*Wopt(:,k) );
%     oSNR(:,k) = abs( Wopt(:,k)'*Rxx(1:M,1:M,k)*Wopt(:,k) ) / abs( Wopt(:,k)'*Rvv(1:M,1:M,k)*Wopt(:,k) );
end
% for k = 1:K/2+1           
%     Rinv = inv( Rvv(1:M,1:M,k) );
%     Wopt(:,k) = Rinv * dvec_s(:,k) / (dvec_s(:,k)' * Rinv * dvec_s(:,k) );
%     oSNR(:,k) = abs( Wopt(:,k)'*Rxx(1:M,1:M,k)*Wopt(:,k) ) / abs( Wopt(:,k)'*Rvv(1:M,1:M,k)*Wopt(:,k) );
% end
%% Senseor Slection
UM = ones(M,M);
matrix = A;
matrix(1:size(A, 1)+1:end) = 1;

aa = zeros(M*K0,1);
u = ones(M,1);
J = K0/Q;
D = diag(sum(A, 2));
L = D-A;
[mm,nn] = eig(L);
[nn1, idx] = sort(mm);
nn2 = nn1(2);
mm2 = mm(:, idx(2));
U0 = repmat(eye(M),1,J);
tic;
t = 0;
ii = 0;
for i = 1:10
    ii = ii+1;
    cvx_solver mosek
    cvx_begin quiet
    cvx_precision high
        variable w(M*K0,1) complex
        variable W0(M,M)
        variable W(M*K0,M) complex       
        expression WQ(M,M,Q)
        sumoNoise = 0;

        for q = 1:Q
            WQ(:,:,q) =U0 * abs(W((q-1)*J*M+1:q*J*M,:));
        end
        
        minimize trace(UM*W0)
        subject to
        for q = 1:Q
            for k = (q-1)*J+1:q*J
                kk = k+1;
                k0 = (k-1)*M +1;           
                w(k0:k0+M-1,1)'*dvec_s(:,kk) == 1;
                W(k0:k0+M-1,:) * dvec_s(:,kk) ==  w(k0:k0+M-1,1);
                W(k0:k0+M-1,:)  == hermitian_semidefinite(M);
%                 real(trace( Rxx(1:M,1:M,kk)*W(k0:k0+M-1,:) )) - real(oSNR(:,kk)*b0*trace( Rvv(1:M,1:M,kk)*W(k0:k0+M-1,:) )  ) >= 0;
                sumoNoise =sumoNoise + real(trace(Rvv(1:M,1:M,kk)*W(k0:k0+M-1,:))) - real(oNoise_min(:,kk)/b0);
            end
            W0>=WQ(:,:,q);
        end
        sumoNoise <= 0;
    cvx_end

    d_w(:,i) = sum(abs(reshape(w, M, K0)),2);


    sel_tem_i2 = sum(abs(reshape(w,M,K0)),2);
    sel_tem_i1 = sum(abs(reshape(aa,M,K0)),2);

    sel_i2 = zeros(M,1); 
    sel_i1 = zeros(M,1); 

%     sel_i2(sel_tem_i2>=1e-3) = 1;
%     sel_i1(sel_tem_i1>=1e-3) = 1;
%     slc = find(sel_i2>=1e-3);
    sel_i2(sel_tem_i2>=1e-1) = 1;
    sel_i1(sel_tem_i1>=1e-1) = 1;
    slc = find(sel_i2>=1e-1);
    P = sum(sel_i2);
    tao = eye(P);
    for j = 1:M-1
        tao = tao+A(slc,slc)^j;
    end
    adjacent_nodes = [];
    
    for ii = 1:length(slc)
        ind = slc(ii);
        
        neighbors = find(matrix(ind, :) > 0);
    
        adjacent_nodes = [adjacent_nodes, neighbors];
    end
    adjacent = [];
    adjacent = unique(adjacent_nodes);
    adjacent = sort(adjacent);
    
    ad = zeros(M,1);
    ad(adjacent)=1;
    
    A1 = A(slc,slc);
    D1 = diag(sum(A1, 2)); 
    L1 = D1-A1;
    [mm,nn] = eig(L1);
    % 对特征值进行排序
    [nn1, idx] = sort(diag(nn));
    if size(idx)== 1
        lamda2 = 1;
    else
        lamda2 = nn1(2);
    end
%     if (lamda2>1e-5) && (sum(abs(sel_tem_i2-sel_tem_i1))<=0.1)
    if (lamda2>1e-5) && (norm(sel_tem_i2-sel_tem_i1)<=1)
%     if sel_i2 == sel_i1
        break;
    end
    aa = w;
    
%     sum(abs(sel_tem_i2-sel_tem_i1))
%     lamda2
    
%     y = zeros(M,1);    
%     for k=1:K0
%         k0 = (k-1)*M +1;           
%         [vec, ev] =  eig(W(k0:k0+M-1,:)) ;
%         [~, ind] = max(diag(ev));
%         y = y +  abs(vec(:,ind)).^2;
%     end
%     y  = 1/ K0 * y;
%     Y = y * y.';
%     UM1 = 1./(Y + 0.05) ;
%     
%     mm2 = mm(:, idx(2));
%     alpha = 5; a0 = 0.05;
    index = find(sel_i2 == 1); %本次选择节点的索引
    
%     u(index) = 0;
    if lamda2>1e-5
        u = ((A*sel_i2).^alpha);
        non_index = setdiff(1:numel(u), index); %非index
        UM0(non_index,non_index) = UM(non_index,non_index);
        UM = 1./(u * u.'+ 1./UM);
        UM(non_index,non_index) = UM0(non_index,non_index);
    else
        u = ((A*sel_i2).^alpha);
        u(index) = 0;
        idx = find(u>0);
        non_idx = setdiff(1:numel(u), idx);
        UM0(non_idx,non_idx)= UM(non_idx,non_idx);
        UM = 1./(u * u.'+ 1./UM);
        UM(non_idx,non_idx) = UM0(non_idx,non_idx);
    end

%     u(index) = u(index)*0.8;
%     u = u.^50;
%     UM = 1./((u * y.') + 1./UM);
%     UU = (u * u.').*(y * y.');
   

end
t = t+toc; 
choose = sum(abs(reshape(w, M, K0)),2);
end