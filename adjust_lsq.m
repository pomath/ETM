function [C, S, So, V, r, dof, cst_pass] = adjust_lsq(Ai,Pi,Li,index, constrains)

    % select the rows that are not zero
    r = ~all(Ai == 0,2);
    
    % don not take into account outliers
    A = Ai(r & index,:);
    L = Li(r & index);
    P = Pi(r & index,r & index);
    
    % to stabilize the inversion
    sc = size(constrains,1);
    A = [A; constrains];
    L = [L; zeros(sc,1)];
    sp = size(P,1);
    for i=1:sc
        P(sp+i,sp+i) = 1;
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
