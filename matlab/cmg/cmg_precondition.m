% function [pfun, H, flag] = cmg_precondition(A,opts)
%
% input:  A: SDD matrix (positive off-diagonals not supported currently)
%
%         opts: optional struct argument with fields   
%         opts.matlab_: if 1 then pfun is slower matlab function, else mex (default mex)
%         opts.valiadate : if 0 then A is assumed to be Laplacian, and
%                          input is not validated (default = 1)
%
%
% output: pfun: preconditioning function for A
% 
%      H: preconditioning hierarchy used internally in pfun
%      flag: failure status variable
%
% For an arbitrary b-side in the null space of A, 
% the system Ax = b can be solved by:
%     x = pcg(A, b, tol,iter, pfun);        t
%
% >> help pcg for documentation
%
%

%% 
% CMG, Copyright (c) 2008-2020  Ioannis Koutis, Gary Miller               
% 
% The CMG solver is distributed under the terms of the GNU General Public  
% Lincense Version 3.0 of the Free Software Foundation.                    
% CMG is also available under other licenses; contact authors for details. 


%% List of functions in this file:
%
% cmg_precondition:    main function
% 
% validate_input
% make_preconditioner
%
% Combinatorial Functions
% steiner_group
% forest_components_
% split_forest_
% level_init
%
% Numerical Functions
% ldl_
% preconditioner_           : matlab prototype
% preconditioner_sd
%
% mex Communication Functions 
% process_hierarchy
% process_sparse_matrix

function [pfun, H, flag] = cmg_precondition(A,opts)

    

    % handle small input 
    if length(A)<500
        pfun = [];
        H = [];
        flag = [];         
        disp('Input matrix is small. Solve Ax=b with A\b');
        return
    end

    global matlab_
    try
        if opts.matlab_ == 1
            matlab_ = 1;
        else
            matlab_ = 0;
        end
    catch
        matlab_ = 0;
    end
    try
        opts.validate;
    catch
        opts.validate = 1;
    end

    %% validate input 

    if opts.validate == 1
        [A_, flag] = validate_input(A);
        if flag == 1
            disp('The input matrix must be symmetric');
            pfun = [];
            return
        elseif flag == 2
            disp('The current version of CMG does not support positive off-diagonals');
            pfun = [];
            return
        end        
    else
        flag = 0;
        A_ = A;
    end
    

    
    S_init = level_init();
    H{1} = S_init;
    
    % A_ is Laplacian
    if length(A_)>length(A)
        H{1}.sd = true;           % original matrix is strongly dominant
    else
        H{1}.sd = false;
    end

    
    %% construct hierarchy 
    loop = true;
    j = 1; 
    h_nnz = 0;
    
    
    while loop
        
        n = length(A_);
      
        %direct method for small size
        if (n<700)            
            H{j}.iterative = 0;
            break
        end
        
                
        dA_ = diag(A_);
        [cI, ~] = steiner_group(A_,dA_);
        nc = max(cI);
        
        H{j}.islast = 0;
        H{j}.A = A_;
        H{j}.invD = 1./(2*dA_);
        H{j}.cI = cI;
        H{j}.nc = nc;
        H{j}.R = sparse(cI, 1:n, ones(n,1), nc, n);  % added for efficiency
        
        % check for full contraction
        if nc == 1
            H{j}.islast = 1;
            H{j}.iterative = 1;
            flag = 1;
            break;
        end

 
        % check for hierarchy stagnation for potentially bad reasons
        h_nnz = h_nnz+nnz(A_);
        if (nc >= (n-1)) || (h_nnz > (5*nnz(H{1}.A))) 
            H{j}.islast = 1;
            H{j}.iterative = 1;
            flag = 3;   % indicates stagnation
            warning('CMG convergence may be slow due to matrix density. Future versions of CMG will eliminate this problem.');
            break;
        end
        
        
        Rt = sparse(double(cI),1:n,1,nc,n);
        A_ = Rt*(H{j}.A)*Rt';
        
        
        j = j+1;
        H{j} = S_init;   % initialize level
    end
    
    % code for last hierarchy level
    if flag == 0                              % no previous error
        H{j}.islast = 1;
        H{j}.iterative = 0;
        H{j}.A = A_(1:n-1,1:n-1); 
        %[L, D, p] = ldl( H{j}.A ,'vector');
        [L, D, p] = ldl_( A_(1:n-1,1:n-1) );
        H{j}.chol.ld = L;
        H{j}.chol.ldT = L';
        H{j}.chol.d = 1./diag(D);
        H{j}.chol.p = p;                        % x = A*y => y(p) = LT\(H{j}.d.*(L\x(p))); 
    end
    
    
    %  determine number of recursive calls between levels
    for k=1:(j-2)
       H{k}.repeat = max(floor(nnz(H{k}.A)/nnz(H{k+1}.A)-1),1);
    end
    if flag == 0
        H{j-1}.repeat = max(floor(nnz(H{j-1}.A)/nnz(H{j}.chol.ld)-1),1);
    else
        H{j-1}.repeat = max(floor(nnz(H{j-1}.A)/nnz(H{j}.A)-1),1);
    end
    
    % H = cell2mat(H);
    
    % create precondition function
    pfun = make_preconditioner(H);

    
end % function

%% make preconditioner
function pfun = make_preconditioner(H)

    global matlab_
    
    if matlab_  && ~H{1}.sd
        pfun = @(b) preconditioner_(H,b,1);  
        return
    end
    
    if ~matlab_  && ~H{1}.sd
        sH = process_hierarchy(H);
        pfun = @(b) mx_preconditioner_(sH,b);
        return
    end
    
    if matlab_  && H{1}.sd
        pfun = @(b) preconditioner_sd(H,b);  
        return
    end
    
    if ~matlab_  && H{1}.sd
        sH = process_hierarchy(H);
        pfun = @(b) preconditioner_sd(sH,b);  
        return
    end
  
end


%% preconditioner_sd
function x = preconditioner_sd(H,b)

    global matlab_

    n = length(b);
    b(n+1) = -sum(b);
    if matlab_
        x = preconditioner_(H,b,1);
    else
        x = mx_preconditioner_(H,b);
    end
    x = x(1:n)-x(n+1); 
end
    
%% initialize hierarchy level
function S = level_init()

    S.sd = true;
    S.islast = 0;
    S.iterative = 1;
    S.A = sparse([]);
    S.invD = 0;
    S.cI = 0;
    S.nc = 0;
    S.R = sparse([]);
    S.repeat = 0;

    % used for last level
    S.chol.ld = sparse([]);
    S.chol.ldT = sparse([]);
    S.chol.d = 0;
    S.chol.p = 0;

end


%% steiner groups
function [cI, nc] = steiner_group(A, dA_)

% input A is Laplacian

    global matlab_

    [M, C1] = min(A);      %C1 represents a unimodal forest
     
    if matlab_
        C = split_forest_(C1);
    else
        C = mx_splitforest_(uint32(C1-1))+1;    %+-1 c-indexing
    end

    efd = abs(M'./dA_);
    if min(efd) < (1/8)   %low effective degree nodes found
        C = update_groups_(A, C, dA_);
    end
    
    [cI, nc, ~] = forest_components_(C);

end


%% forest components 
% connected components in unimodal forest
function [cI, nc, csizes] = forest_components_(C)
% input C represents a unimodal tree


    n = length(C);
    cI = zeros(n,1);
    csizes = zeros(n,1);
    buffer = zeros(100,1);   

    ccI = 1;
    for j=1:n
        bufferI = 1;
        jwalk = j;
        % tree walk
        while cI(jwalk)==0
           
            cI(jwalk) = ccI;
            buffer(bufferI) = jwalk; 
            bufferI = bufferI+1;
            
            if bufferI == length(buffer)     % for memory efficiency
                buffer_ = zeros(min(2*length(buffer),n),1); 
                buffer_(1:length(buffer))  = buffer;
                buffer = buffer_;
                clear buffer_;
            end
                
            jwalk = C(jwalk);
        end %while
        
        bufferI = bufferI-1;
        en = C(jwalk);  % end node 
        if (cI(en) ~= ccI)
            cI(buffer(1:bufferI)) = cI(en);          
        else
            ccI = ccI+1;
        end
        csizes(en) = csizes(en) + bufferI;

    end %for
    if csizes(ccI) == 0
        ccI = ccI-1;
    end
    nc = ccI;
    csizes = csizes(1:ccI);

end



%% update groups based on nodes with low effective degree
function C1 = update_groups_(A, C, dA_)

% input A is Laplacian

    n = length(C);
    B = zeros(n,1);

    %B(j) is the total tree weight incident to node j
    for j=1:n    
        k = C(j);
        if k ~=j
            a = A(j,k);
            B(j)=  a+B(j);
            B(k) = a+B(k);
        end
    end

    ndx  = find((B./dA_)>-0.125);  % negative sign because A is laplacian 
    C(ndx) = uint32(ndx); 
    
    C1 = C;
end

%% input validation
function [A0, flag] = validate_input(A)

% A0 is Laplacian, possibly extended by one coordinate when A is SD
% failure_flag code indicates issue with input validity

    flag = 0;

    % check symmetry 
    if not(issymmetric(A))
        flag = 1; 
    end

    % detect strict dominance
    n = length(A);
    dA = diag(A);
    sA = sum(A)'; 
    sAp = (sA+abs(sA))/2;
    sd = sAp./dA>1e-13;        % if sd(i) =1 then row i is numerically strictly dominant


    [i,j,v] = find(A);

    % check for positive off-diagonals that are not currently supported    
    if max( v(not(i==j)) ) >0
        flag = 2;
    end
    
    % augment by extra coordinate if strictly dominant    
    if max(sd)
       ex_v = -sAp(sd);  
       exd = length(ex_v);
       exd_c = find(sd);       %get coordinates
       
       i = [i; (n+1)*ones(exd,1); exd_c; n+1]; 
       j = [j; exd_c; (n+1)*ones(exd,1); n+1];
       v = [v; ex_v; ex_v; -sum(ex_v)];  
       A0 = sparse(i,j,v,n+1,n+1);
    else
       A0 = A;
    end
    

end

%% decompose unimodal forest into low conductance components
function C = split_forest_(C1)

    n = length(C1);
    
    C = C1;
    ancestors = zeros(n,1);
    indegree = zeros(n+2,1); 
    visited = zeros(n,1,'logical');
  
    
    walkbuffer = zeros(20,1);
    newancestorbuff = zeros(20,1);
    
    % compute indegrees
    for j=1:n 
        indegree(C(j)) = indegree(C(j))+1; 
    end
        
    % partition into clusters of small diameter
    for j=1:n
        jwalk = j;
        startwalk = true;
        
        while (startwalk && (indegree(jwalk) == 0) && ~visited(jwalk)  )
            startwalk = false;
            ancestors_in_path = 0;   % change over C-CMG !!
            k=1;
            walkbuffer(k) = jwalk;
            newancestorbuff(k)=0;
            while (k<=6 || visited(jwalk))    %+1 for c-indexing adjust
                jwalk = C(jwalk);
                walkterminated = (jwalk == walkbuffer(k)) || ((k>1) && (jwalk == walkbuffer(k-1)));
                if (walkterminated)
                    break; %while
                end
                k = k+1;
                walkbuffer(k) = jwalk;
                if visited(jwalk)
                    newancestorbuff(k) = ancestors_in_path;
                else 
                    ancestors_in_path = ancestors_in_path+1;
                    newancestorbuff(k) = ancestors_in_path;   
                end
            end
                
            if (k>6) % large diameter - cut
                middlek = ceil(k/2);
                C(walkbuffer(middlek)) = walkbuffer(middlek);  % cut middle edge
                indegree(walkbuffer(middlek+1)) = indegree(walkbuffer(middlek+1))-1; % update indegree 
                for ik=(middlek+1):k 
                    ancestors(walkbuffer(ik)) = ancestors(walkbuffer(ik))-ancestors(walkbuffer(middlek));
                end
                for ik=1:middlek    % update ancestors and visited flag   ??
                    visited(walkbuffer(ik)) = true;
                    ancestors(walkbuffer(ik)) = ancestors(walkbuffer(ik))+newancestorbuff(ik);
                end
                jwalk = walkbuffer(middlek+1);  %set first vertex in new walk 
                startwalk = true; 
            end %end cut procedure
            
            if ~startwalk   % commit walk changes          
                for ik=1:k
                    ancestors(walkbuffer(ik)) = ancestors(walkbuffer(ik))+newancestorbuff(ik);                    
                    visited(walkbuffer(ik))=true;
                end
            end
        end   % outer while
    end 
    
    
        
    % tree partition into clusters of high conductance 
    for j=1:n
        jwalk = j;
        startwalk = true;
        
        while (startwalk && (indegree(jwalk) == 0))
            startwalk = false;            
            jwalkb = jwalk;
            cut_mode = false;
            
            while true
                jwalka = C(jwalk);
                walkterminated= ((jwalka == jwalk) || (jwalka == jwalkb));
                if (walkterminated)
                    break; %while
                end
                    

                if (~cut_mode && (ancestors(jwalk) >2) && ((ancestors(jwalka)-ancestors(jwalk))>2)) % possibly low conductance - make cut
                    C(jwalk) = jwalk; % cut edge 
                    indegree(jwalka) = indegree(jwalka)-1;
                    removed_ancestors = ancestors(jwalk);
                    new_front = jwalka;
                    cut_mode = true;
                end  % end making cut
                
                jwalkb = jwalk;
                jwalk = jwalka;
                if cut_mode
                    ancestors(jwalk) = ancestors(jwalk) - removed_ancestors;
                end
            end      
            if cut_mode
                startwalk = true;
                jwalk = new_front;
            end
        end 
    end 
            
 
end % split_forest_
            

%% ldl factorization with min degree heuristic
function [L,D,p] = ldl_(A)

    n = length(A);
    D = zeros(n,1);

    % trivial min degree heuristic: sorting degrees
    [iA, jA] = find(A);
    Ao = sparse(iA,jA,ones(length(iA),1),n,n);
    dA = sum(Ao);
    [~, p] = sort(dA);
    
    A = A(p,p);
    L = speye(n);
    for I=1:n
        vi=A(I+1:n,I); di=A(I,I); D(I)=di;
        L(I+1:n,I) = vi/di;
        Z = vi*vi'/di;
        [Iz, Jz, V] = find(Z); 
        Z = sparse(Iz+I,Jz+I,V,n,n);
        A = A -Z;
    end
    
    D = sparse(1:n,1:n,D, n,n);

end
            


%% preconditioner
function x = preconditioner_(H,b,level)

    n = length(H{level}.A);

    % base case I: last level direct
    if (H{level}.islast && ~H{level}.iterative)

        b = b(1:n);
        x = zeros(n+1,1);

        L = H{level}.chol.ld; 
        LT = H{level}.chol.ldT;
        di = H{level}.chol.d;
        p = H{level}.chol.p;

        x(p) = LT\(di.*(L\b(p)));  % back substitution needed here
        
        return
    end
    
    % base case II: last level iterative
    if (H{level}.islast && H{level}.iterative)
        invD = H{level}.invD;
        x = invD.*b;
        return
    end

    % main cycle
    
        % unpack variables for efficiency
        invD = H{level}.invD;
        A = H{level}.A;
        repeat = H{level}.repeat;
        cI = H{level}.cI;
        R = H{level}.R;

    
    for j=1:repeat
                   
        % Jacobi pre-smooth 
        if j==1 
            x = invD.*b;            % previous x = 0
        else
            x = x + invD.*(b-A*x);
        end
        
        
        % residual 
        r = b - A*x; 

        % projection of residual onto smaller level
        % b_small = zeros(nc,1);
        % for k=1:n
        %     b_small( cI(k) ) = b_small( cI(k) )+ r(k);
        % end
        b_small = R*r;

        z = preconditioner_(H,b_small,level+1);
        
        %interpolation
        z = z(cI); 
        
        x = x+z;
        
        % % Jacobi post-smooth 
        x = x + invD.*(b-A*x);
    
    end
    
    
end


%% process_sparse_matrix
function sa = process_sparse_matrix(A,sym,pre)

    if not(or(sym==1,sym==0))
        error('first input must be 1 or 0');
    end

    if not(or(pre==1,pre==2))
        error('second input must be 1 or 2');
    end


    n = length(A);
    if (sym == 0)
        A=A';
    end
    if (sym == 1)
        A = tril(A);
    end    


    [i, j, v] = find(A);

    [~, col_j] = unique(j,'first');
    col_j = uint32(col_j-1);
    col_j(n+1)= uint32(length(v));
    row_i = uint32(i-1);


    sa.a = v;
    sa.ia = col_j; 
    sa.ja = row_i;
    sa.issym = sym;
    sa.n = n;

end

%% process hierarchy
function H = process_hierarchy(H)


    nH = length(H)-1;

    for j=1:(nH)
        H{j}.A = process_sparse_matrix(H{j}.A,1,2);

        % uint32 and prepare for C-indexing
        H{j}.cI = uint32(H{j}.cI-1);


        % auxiliary work space for preconditioner
        if  or(j<nH, H{j}.iterative==1)
            H{j}.lws1 = zeros(H{j}.A.n+1,1);
            H{j}.lws2 = zeros(H{j}.A.n+1,1);
            H{j}.sws1 = zeros(H{j}.nc+1,1);
            H{j}.sws2 = zeros(H{j}.nc+1,1);
            H{j}.sws3 = zeros(H{j}.nc+1,1);        
        end
    end
    
    % last level
    H{end}.A = process_sparse_matrix(H{end}.A,1,2);
    
    if ~H{end}.iterative
        
        ld = H{end}.chol.ld;
        d = 1./H{end}.chol.d;            % re-inverting
        ld = ld-eye(length(d))+diag(d);  %d is stored in diagonal
        p = H{end}.chol.p;
        invp(p) = 1:length(p);

        H{end}.chol.ld = process_sparse_matrix(ld,0,2);
        H{end}.chol.ldT = process_sparse_matrix(ld',0,2);
        H{end}.chol.p = uint32(p-1);
        H{end}.chol.invp = uint32(invp-1);
    end

end
