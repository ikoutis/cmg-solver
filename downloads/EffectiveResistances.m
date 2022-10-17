function [eff_res] = EffectiveResistances(elist,e,w,tol,epsilon,type_,pfun_)
% function [eff_res] = EffectiveResistances(elist,w,tol,epsilon,type_,pfun_)
%
% calculate effective resistances in a given graph:
% the graph is viewed as an electrical resistive network
% where the edge weights correspond to capacitances.
%
%
% Input:
%   elist: list of edges whose effective resistance is needed [K x 2]
%   e: edges of the graph [M x 2]
%   w: weights in the graph [M x 1]
%   [tol]: specify the relative error in the linear system solution
%            1e-4 <default>
%   [epsilon]: controls accuracy in the computation of effective resistance
%              decreasing epsilon increases accuracy and computation time
%              1 <default>
%   [type_]: optional string
%           'slm' <default>: this version use less memory.
%           'spl': this implements the Spielman-Srivastava algorithm
%           'org': find the exact effective resistances. Works only for
%                  small networks
%   [pfun_]: cmg_sdd(L), if already computed.
% Output:
%   eff_res: vector of effective resistances for the edges in elist
%
%
% Example: The path graph.
%   E = [(1:49)' (2:50)']; w = ones(length(E),1); elist = [1 50];
%   er = EffectiveResistances(elist,E,w,1e-5,1,'org')
%
% Richard Garcia-Lebron
%
    %% input Validation
    [m,two] = size(e);
    [~,two1] = size(elist);
    if nargin > 7
        error('Too many arguments, see help EffectiveResistancesJL.');
    elseif nargin < 3
        error('More arguments needed, see help EffectiveResistancesJL.');
    elseif nargin == 5
        type_ = 'slm';
    elseif nargin == 4
        type_ = 'slm';
        tol = 10^-4; %tolerance for the cmg solver
        epsilon = 1;
    elseif nargin == 3
        type_ = 'slm';
        tol = 10^-4; %tolerance for the cmg solver
        epsilon = 1;
    elseif two ~= 2 ||  two1 ~= 2
        estring = ['The first or the second input have wrong' ...
                    ' column dimension should be [M  x 2].'];
        error(estring);
    end

    %% Creating data for effective resitances
    numIterations = 300; %iteration for the cmg solver
    tolProb = 0.5;%use to create the johnson lindestraus projection matrix
    n = max(max(e));%number of node in the graph
    A = sparse([e(:,1);e(:,2)],[e(:,2);e(:,1)],[w;w],n,n); %adjacency matrix
    L = diag(sum(abs(A),2)) - A; %laplacian matrix
    clear 'A' %adjacensy matrix no needed
    B = sparse([1:m 1:m],[e(:,1) e(:,2)],[ones(m,1) -1*ones(m,1)]);
    [m,n] = size(B);
    W = sparse(1:length(w),1:length(w),w.^(1/2),length(w),length(w));
    scale = ceil(log2(n))/epsilon;
    %% Finding the effective resitances
    if strcmp(type_,'slm')
        if nargin == 7 % Creating the preconditioned function
            pfun = pfun_;
        elseif (length(L) > 600)
            pfun = cmg_sdd(L);
        end
        [rows,~,~] = size(elist);
        eff_res = zeros(1,rows);
        optimset('display','off');
        if (length(L) > 600) % bigger graphs
            for j=1:scale
                ons = (rand(1,m) > tolProb);
                ons = ons - not(ons);
                ons = ons./sqrt(scale);
                [Z flag]= pcg(L,(ons*(W)*B)',tol,numIterations,pfun);
                if flag > 0
                    error(['PCG FLAG: ' num2str(flag)])
                end
                Z = Z';
                eff_res = eff_res + (abs((Z(elist(:,1))-Z(elist(:,2)))).^2);
            end
        else % smaller graphs
            for j=1:scale
                ons = (rand(1,m) > tolProb);
                ons = ons - not(ons);
                ons = ons./sqrt(scale);
                [Z flag]= pcg(L,(ons*(W)*B)',tol,numIterations);
                if flag > 0
                    error(['PCG FLAG: ' num2str(flag)])
                end
                Z = Z';
                eff_res = eff_res + (abs((Z(elist(:,1))-Z(elist(:,2)))).^2); % the second "loop"
            end
        end
        eff_res = eff_res';
    elseif strcmp(type_,'spl')
        Q = sparse((rand(scale,m)) > tolProb);
        Q = Q - not(Q);
        Q = Q./sqrt(scale);
        SYS = sparse(Q*(W)*B); % Creation the system
        if nargin == 7 % Creating the preconditioned function
            pfun = pfun_;
        elseif (length(L) > 600)
            pfun = cmg_sdd(L);
        end
        optimset('display','off');
        %% Solving the systems
        if (length(L) > 600) % bigger graphs
            for j=1:scale
                [Z(j,:) flag] = pcg(L,SYS(j,:)',tol,numIterations,pfun);
                if flag > 0
                    error(['PCG FLAG: ' num2str(flag)])
                end
            end
        else % smaller graphs
            for j=1:scale
                [Z(j,:) flag] = pcg(L,SYS(j,:)',tol,numIterations);
                if flag > 0
                    error(['PCG FLAG: ' num2str(flag)])
                end
            end
        end
        eff_res = sum(((Z(:,elist(:,1))-Z(:,elist(:,2))).^2),1)';
    elseif strcmp(type_,'org')
        [m,~] = size(elist);
        B = sparse([1:m 1:m],[elist(:,1) elist(:,2)],[ones(m,1) -1*ones(m,1)],m,length(L));
            if length(L) > 600 % bigger graphs
                if nargin == 7
                    pfun = pfun_;
                else
                    pfun = cmg_sdd(L);
                end
                optimset('display','off');
                for j=1:m
                    [Z flag]= pcg(L,B(j,:)',1e-10,numIterations,pfun);
                    if flag > 0
                        error(['PCG FLAG: ' num2str(flag)])
                    end
                    eff_res(j) = B(j,:)*Z;
                end
            else % smaller graphs
                opts.type = 'nofill'; opts.michol = 'on';
                for j=1:m
                    [Z,flag] = pcg(L,B(j,:)',1e-10,numIterations);
                    eff_res(j) = B(j,:)*Z;
                end
            end
        eff_res = eff_res';
    end
end

