function varargout = CHPSO_func(fhd,numIter,param,varargin)

format short e

% make copies of parameters
NP  = param.ParicleNumber;
acc = param.AccelerationConstants;
iwt = param.InteriaWeight;
dim = param.Dimension;
Ub  = param.Ub;
Lb  = param.Lb;


mv = repmat(0.5*(Ub-Lb),1,dim);

w_max = iwt(1);
w_min = iwt(2);

pos_a = Lb+(Ub-Lb).*rand(NP,dim);
pos_b = Lb+(Ub-Lb).*rand(NP,dim);
pos_c = Lb+(Ub-Lb).*rand(NP,dim);
pos_d = Lb+(Ub-Lb).*rand(NP,dim);

vel_a = Lb+(Ub-Lb).*rand(NP,dim);
vel_b = Lb+(Ub-Lb).*rand(NP,dim);
vel_c = Lb+(Ub-Lb).*rand(NP,dim);
vel_d = Lb+(Ub-Lb).*rand(NP,dim);

for j = 1:NP
    val_a(j,1) = feval(fhd,pos_a(j,:),varargin{:});
    val_b(j,1) = feval(fhd,pos_b(j,:),varargin{:});
    val_c(j,1) = feval(fhd,pos_c(j,:),varargin{:});
    val_d(j,1) = feval(fhd,pos_d(j,:),varargin{:});
end
localPos_a = pos_a; localVal_a = val_a;
localPos_b = pos_b; localVal_b = val_b;
localPos_c = pos_c; localVal_c = val_c;
localPos_d = pos_d; localVal_d = val_d;

[gval_a, idx_a] = min(val_a(:,1));  globalPos_a = pos_a(idx_a,:);
[gval_b, idx_b] = min(val_b(:,1));  globalPos_b = pos_b(idx_b,:);
[gval_c, idx_c] = min(val_c(:,1));  globalPos_c = pos_c(idx_c,:);
[gval_d, idx_d] = min(val_d(:,1));  globalPos_d = pos_d(idx_d,:);

G = [globalPos_a; globalPos_b; globalPos_c; globalPos_d];
[gVal, idx] = min([gval_a; gval_b; gval_c; gval_d]);
globalPos = G(idx,:);

%   best record
BestChart = [];
BestChart = [BestChart; gVal];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rfc = 0; RFC = 2;
for i = 2:numIter
    fit_avg = (sum(val_a(:,1))+sum(val_b(:,1))+sum(val_c(:,1))+sum(val_d(:,1)))/(4*NP);
    fit_m = [min(val_a(:,1)),min(val_b(:,1)),min(val_c(:,1)),min(val_d(:,1))];
    fit_min = min(fit_m);
    for j = 1:NP
        if (val_a(j,1)+val_b(j,1))/2>=fit_avg
            w = w_max;            
        elseif (val_a(j,1)+val_b(j,1))/2>fit_avg/2&&(val_a(j,1)+val_b(j,1))/2<fit_avg            
            w = w_min+(w_max-w_min)*((val_a(j,1)+val_b(j,1))/2-fit_avg/2)/(fit_avg/2-fit_min);
        else
            w = w_max-(w_max-w_min)*i/numIter;            
        end
        
        vel_a(j,:) =   w.*vel_a(j,:)+acc(1).*rand.*(localPos_a(j,1:dim)-pos_a(j,1:dim))...
                    +acc(2).*rand.*(globalPos(1:dim)-pos_a(j,1:dim));
        vel_a(j,:) =  (vel_a(j,:)>mv).*mv+(vel_a(j,:)<=mv).*vel_a(j,:);       
        vel_a(j,:) =  (vel_a(j,:)<(-mv)).*(-mv)+(vel_a(j,:)>=(-mv)).*vel_a(j,:);
        
        pos_a(j,1:dim) = pos_a(j,1:dim)+vel_a(j,1:dim);     
        val_a(j,1) = feval(fhd,pos_a(j,1:dim),varargin{:});
        
        vel_b(j,:) =   w.*vel_b(j,:)+acc(1).*rand.*(localPos_b(j,1:dim)-pos_b(j,1:dim))...
                    +acc(2).*rand.*(globalPos(1:dim)-pos_b(j,1:dim));
        vel_b(j,:) =  (vel_b(j,:)>mv).*mv+(vel_b(j,:)<=mv).*vel_b(j,:);       
        vel_b(j,:) =  (vel_b(j,:)<(-mv)).*(-mv)+(vel_b(j,:)>=(-mv)).*vel_b(j,:);
        
        pos_b(j,1:dim) = pos_b(j,1:dim)+vel_b(j,1:dim);     
        val_b(j,1) = feval(fhd,pos_b(j,1:dim),varargin{:});
        
        gama1 = val_a(j,1);
        gama2 = val_b(j,1);
        gama = val_a(j,1)+val_b(j,1);
        
        vel_c(j,:) =   w.*((gama.*vel_a(j,:)/gama1+gama.*vel_b(j,:)/gama2)+vel_c(j,:))...
                    +acc(1).*rand.*(localPos_c(j,1:dim)-pos_c(j,1:dim))...
                    +acc(2).*rand.*(globalPos(1:dim)-pos_c(j,1:dim));
        vel_c(j,:) =  (vel_c(j,:)>mv).*mv+(vel_c(j,:)<=mv).*vel_c(j,:);       
        vel_c(j,:) =  (vel_c(j,:)<(-mv)).*(-mv)+(vel_c(j,:)>=(-mv)).*vel_c(j,:);
        
        pos_c(j,1:dim)  = pos_c(j,1:dim)+vel_c(j,1:dim);
        val_c(j,1)    = feval(fhd,pos_c(j,1:dim),varargin{:});
    
        vel_d(j,:)=(vel_a(j,:)+vel_b(j,:)-vel_c(j,:));
        vel_d(j,:) =  (vel_d(j,:)>mv).*mv+(vel_d(j,:)<=mv).*vel_d(j,:);       
        vel_d(j,:) =  (vel_d(j,:)<(-mv)).*(-mv)+(vel_d(j,:)>=(-mv)).*vel_d(j,:);
        pos_d(j,1:dim) = pos_d(j,1:dim)./6+localPos_d(j,1:dim)./3+globalPos(1,1:dim)./2+vel_d(j,1:dim);
        val_d(j,1) = feval(fhd,pos_d(j,1:dim),varargin{:});

        if (val_a(j,1)<localVal_a(j,1))
            localPos_a(j,:) = pos_a(j,:);
            localVal_a(j,1) = val_a(j,1);
        end

        if (val_b(j,1)<localVal_b(j,1))
            localPos_b(j,:) = pos_b(j,:);
            localVal_b(j,1) = val_b(j,1);
        end

        if (val_c(j,1)<localVal_c(j,1))
            localPos_c(j,:) = pos_c(j,:);
            localVal_c(j,1) = val_c(j,1);
        end

        if (val_d(j,1)<localVal_d(j,1))
            localPos_d(j,:) = pos_d(j,:);
            localVal_d(j,1) = val_d(j,1);
        end
    end

    [gval_a, idx_a] = min(val_a(:,1));
    [gval_b, idx_b] = min(val_b(:,1));
    [gval_c, idx_c] = min(val_c(:,1));
    [gval_d, idx_d] = min(val_d(:,1));
    globalPos_a = pos_a(idx_a,:);
    globalPos_b = pos_b(idx_b,:);
    globalPos_c = pos_c(idx_c,:);
    globalPos_d = pos_d(idx_d,:);
    
    [gVal_tmp idx] = min([gval_a; gval_b; gval_c; gval_d]);   
    G = [globalPos_a; globalPos_b; globalPos_c; globalPos_d];

    if gVal_tmp<gVal
        globalPos = G(idx,:);
        gVal = gVal_tmp;
    end
    
    BestChart = [BestChart; gVal];

    rfc = rfc+1;
    if rfc >= RFC
        rfc = 0;
        [wval_a, b_index1] = max(val_a(:,1));
        [wval_b, b_index2] = max(val_b(:,1));
        [wval_c, b_index3] = max(val_c(:,1));
        [wval_d, b_index4] = max(val_d(:,1));
        pos_a(b_index1,:) = -globalPos;
        pos_b(b_index2,:) = -globalPos;
        pos_c(b_index3,:) = -globalPos;
        pos_d(b_index4,:) = -globalPos;        
    end    
end

varargout{1} = BestChart;




