function [C, S, So, V, r, dof, cst_pass, sigma, index] = adjust_lsq(Ai,Pi,Li,constrains_h,constrains_s)

    limit = 2.5;
    
    % select the rows that are not zero
    r = ~all(Ai == 0,2);
    
    % don not take into account outliers
    A  = Ai(r,:);
    Ai = Ai(r,:);
    
    L  = Li(r);
    Li = Li(r);
    
    P  = Pi(r,r);
    Pi = Pi(r,r);
    
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
    
    % to stabilize the inversion, add the conditions equations
    sc = size(constrains,1);
    A = [A; constrains];
    L = [L; zeros(sc,1)];
    
    cst_pass = false;
    iter = 0;
    factor = 1;
    
    while ~cst_pass & iter <= 10
        % REBUILD the P matrix. Inside the loop to reweigh after adjustment
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

        dof = (size(Ai,1)-size(Ai,2));

        So = sqrt(V'*Pi*V/dof);

        x = So.^2.*dof;

        % find the a priori sigma for the observations
        factor = factor.*So;
        % find how many sigmas away the outliers are
        s = V./factor;
            
        % careful! This function returns the opposite value of alpha as on
        % Leick, page 143
        X1 = chi2inv(1-0.05/2,dof);
        X2 = chi2inv(0.05/2,dof);

        if x < X2 || x > X1
            % if it falls in here it's because it didn't pass the Chi2 test
            cst_pass = false;
            
            if So < 1
                % weights are too pesimistic, just inform the user
                fprintf(char(hex2dec('25B2')));
            else
                % weights are too optimistic, just inform the user
                fprintf(char(hex2dec('25BC')));
            end
            
            % reweigh by Mike's method of equal weight until 2 sigma
            f = ones(size(V));
            f(s > limit) = 1./(10.^(limit - s(s > limit)));
            % do not allow sigmas > 100 m, which is basicaly not putting
            % the observation in. Otherwise, due to a model problem
            % (missing jump, etc) you end up with very unstable inversions
            f(f > 100) = 100;
            
            P = diag(1./((factor.*f).^2));
            Pi = P;
        else
            cst_pass = true;
        end
        
        iter = iter + 1;
    end
    
    %%%%%%%%%%%% statistics %%%%%%%%%
    S = inv(Ai'*Pi*Ai);
    sigma = diag(1./sqrt(Pi));
    
    % mark observations with sigma > limit
    index = true(size(V));
    index(s > limit) = 0;
    
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
