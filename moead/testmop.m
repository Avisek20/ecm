function mop = testmop( testname, dimension, domain, K, data, ecmMode )
%Get test multi-objective problems from a given name.
%   The method get testing or benchmark problems for Multi-Objective
%   Optimization. The implemented problems included ZDT, OKA, KNO.
%   User get the corresponding test problem, which is an instance of class
%   mop, by passing the problem name and optional dimension parameters.



mop=struct('name',[],'od',[],'pd',[],'domain',[],'func',[]);
switch lower(testname)
    case 'kno1'
        mop=kno1(mop);
    case 'zdt1'
        mop=zdt1(mop, dimension);
    case 'ecm'
        mop = ecm(mop, dimension, domain, K, data, ecmMode);
    otherwise
        error('Undefined test problem name');
end
end

% Entropy kmeans clustering
function p = ecm(p, dim, domain, K, data, ecmMode)
p.name='ECM';
p.pd = dim;
p.domain= domain;
datasetDim = dim / K;
global sigma;
if ecmMode == 1
    p.func = @(x) entropy_kmeans(x, datasetDim, K, data);
    p.od = 2;
else
    p.func = @(x) entropy_kmeans_indiv(x, datasetDim, K, data);
    p.od = K + 1;
end

    function f = entropy_kmeans_indiv(x, datasetDim, K, data)
        
        z = getCentroidFromChromosome(x, datasetDim); % Obtain centroids
        N = size(data, 1);
        CE = 0; CK = zeros(K, 1);
        u = zeros(N, K);
        %beta = 10;
        
        for i=1:N
            
            for j = 1:K
				u(i,j) = exp( -(1/(2*(sigma^2))) * (sum((data(i,:)-z(j,:)).^2)));
			end
			tmp_sum = sum(u(i,:));
			if tmp_sum == 0
				u(i,:) = 1/K;
			else
				u(i,:) = u(i,:)/sum(u(i,:));
			end
   
			% Summing up the objective function values
			for j = 1:K
				CK = CK + u(i,j)*(sum((data(i,:)-z(j,:)).^2));
				if u(i,j) > 1e-16
					CE = CE + u(i,j)*log(u(i,j))/log(2);
				end
			end
        end
        
        f = [CK; CE];
    end

    function f = entropy_kmeans(x, datasetDim, K, data)
        
        z = getCentroidFromChromosome(x, datasetDim); % Obtain centroids
        N = size(data, 1);
        CE = 0; CK = 0;
        u = zeros(N, K);
        %beta = 10;
        
        for i=1:N
            
            for j = 1:K
				u(i,j) = exp( -(1/(2*(sigma^2))) * (sum((data(i,:)-z(j,:)).^2)));
			end
			tmp_sum = sum(u(i,:));
			if tmp_sum == 0
				u(i,:) = 1/K;
			else
				u(i,:) = u(i,:)/sum(u(i,:));
			end
   
			% Summing up the objective function values
			for j = 1:K
				CK = CK + u(i,j)*(sum((data(i,:)-z(j,:)).^2));
				if u(i,j) > 1e-16
					CE = CE + u(i,j)*log(u(i,j))/log(2);
				end
			end
        end
        
        f = [CK; CE];
    end
end

%KNO1 function generator
function p=kno1(p)
p.name='KNO1';
p.od = 2;
p.pd = 2;
p.domain= [0 3;0 3];
p.func = @evaluate;

%KNO1 evaluation function.
    function y = evaluate(x)
        y=zeros(2,1);
        c = x(1)+x(2);
        f = 9-(3*sin(2.5*c^0.5) + 3*sin(4*c) + 5 *sin(2*c+2));
        g = (pi/2.0)*(x(1)-x(2)+3.0)/6.0;
        y(1)= 20-(f*cos(g));
        y(2)= 20-(f*sin(g));
    end
end

%ZDT1 function generator
function p=zdt1(p,dim)
p.name='ZDT1';
p.pd=dim;
p.od=2;
p.domain=[zeros(dim,1) ones(dim,1)];
p.func=@evaluate;

%KNO1 evaluation function.
    function y=evaluate(x)
        y=zeros(2,1);
        y(1) = x(1);
        su = sum(x)-x(1);
        g = 1 + 9 * su / (dim - 1);
        y(2) =g*(1 - sqrt(x(1) / g));
    end
end


