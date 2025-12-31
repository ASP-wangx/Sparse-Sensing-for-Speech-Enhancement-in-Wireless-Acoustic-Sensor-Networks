function W0  = EEW_opt_symm(A)


    M = size(A,2);
    A  = A + eye(M);
    %cvx_solver SeDuMi
    cvx_begin quiet
    cvx_precision high
    
        variable W(M,M) symmetric
        variable s
        
        Z1 =  W - 1/M*ones(M,1)*ones(1,M) + s*eye(M);
        Z2 = s*eye(M) -  W + 1/M*ones(M,1)*ones(1,M) ;
        
        minimize s
%         minimize s
        subject to
%             s  <= 1;
            W*ones(M,1) == ones(M,1);
            Z1 == hermitian_semidefinite(M);
            Z2 == hermitian_semidefinite(M);
            % 
            for i=1:M
                for j=1:M 
                    if A(i,j)==0
                    W(i,j) == 0;
                    end
                end
            end

    cvx_end
    
    if isinf(cvx_optval)
        disp('cvx result is inf');
        return;
    end
    
    W0 = full(W);
    
end