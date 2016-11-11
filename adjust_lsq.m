function [C, S, So, V, r, dof, cst_pass] = adjust_lsq(Ai,Pi,Li,index, constrains_h,constrains_s)

    % select the rows that are not zero
    r = ~all(Ai == 0,2);
    
    % don not take into account outliers
    A = Ai(r & index,:);
    L = Li(r & index);
    P = Pi(r & index,r & index);
    
    % "Smart" stabilization inversion
    % try 4 combinations of A:
    % 1) without any condition equations
    % 2) with jumps condition equations
    % 3) with sine and coside condition equations
    % 4) with jump and sine and cosine condition equations
    % pick the design matrix that yields the best stability
    if ~isempty(constrains_h) | ~isempty(constrains_s)
        constrains = find_stable_A(A,constrains_h,constrains_s);
    else
        constrains = [];
    end
    % to stabilize the inversion
    sc = size(constrains,1);
    A = [A; constrains];
    L = [L; zeros(sc,1)];
    sp = size(P,1);
    if sc > 0
        for i=1:sc
            P(sp+i,sp+i) = 1;
        end
    end
    % invert for the parameters
    C = (A'*P*A)\A'*P*L;

    % use the input A to get the information about outliers too
    V = Li - Ai*C;
    %%%%%%%%%%%% statistics %%%%%%%%%
    S = inv(A'*P*A);

    dof = (size(A,1)-size(A,2));
    
    So = sqrt(V'*Pi*V/dof);
    
    x = So.^2.*dof;
    
    % careful! This function returns the opposite value of alpha as on
    % Leick, page 143
    X1 = chi2inv(1-0.05/2,dof);
    X2 = chi2inv(0.05/2,dof);
    
    if x < X2 || x > X1
        % if it falls in here it's because it didn't pass the Chi2 test
        cst_pass = false;
    else
        cst_pass = true;
    end

end

function rcondeq = find_stable_A(A, cond_h, cond_s)
    
    % build the different combinations of
    % condition equations
    condeq{1} = cond_h;
    condeq{2} = cond_s;
    condeq{3} = [cond_s; cond_h];
    
    condnum(1) = cond(A);
    
    for i = 1:length(condeq) 
        condnum(i+1) = cond([A; condeq{i}]);
    end
    
    % find the minimum
    [~,i] = min(condnum);
    if i == 1
        rcondeq = [];
    else
        rcondeq = condeq{i-1};
    end
end
