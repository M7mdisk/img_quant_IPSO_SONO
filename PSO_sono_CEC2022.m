function [gbestval,ccurve, dcurve,gbest] = PSO_sono_CEC2022(ps, nfe_max, Xmin, Xmax, D,fhd,fid)

    gen_max = nfe_max/(ps+1);
    w_max = 0.9;  w_min = 0.4;
    c3 = D/ps*0.01;
    neighbor(1,:)=[ps,2];
    for i=2:ps-1
        neighbor(i,:)=[i-1,i+1];
    end
    neighbor(ps,:)=[ps-1,1];
    pos = Xmin + (Xmax-Xmin).*rand(ps,D);
    eval = feval(fhd,pos',fid);
    nfe = ps;
    pbest = pos;
    pbestval = eval;
    [gbestval,gbestid] = min(pbestval);
    gbest = pbest(gbestid,:);
    gbestrep = repmat(gbest,ps,1);
    Vmax = 30;
    Vmin = -Vmax;
    vel = Vmin+(Vmax-Vmin).*rand(ps,D);

    iter = 2;
    ratio = 0.5;

    %% curves
    ccurve = zeros(1,nfe_max);
    ccurve(1:nfe) = min(pbestval);
    %diversity
    dcurve = zeros(1,nfe_max);
    % Determine the initial diversity of the archive
    prum = mean(pos);
    prum_mat=repmat(prum,ps,1);
    diam2 = (sum( sum( (pos-prum_mat).*(pos-prum_mat) ) ) / ps);
    diversity = sqrt(diam2);
    dcurve(1:nfe) = diversity;
    %%

    while nfe<nfe_max

        %% curves
        old_nfe = nfe;
        %%
        w = w_max - (w_max-w_min)*iter/gen_max;
        [~,index]=sort(eval,'descend');
        pos = pos(index,:);
        vel = vel(index,:);
        pbest = pbest(index,:);
        pbestval = pbestval(index);


        center=ones(ps,1)*mean(pos);
        winidxmask = repmat((1:ps)',1,D);
        winidx = winidxmask + ceil(rand(ps,D).*(ps-winidxmask));
        pwin = pos;
        for j = 1:D
            pwin(:,j) = pos(winidx(:,j),j);
        end


        c1 =(2.5-2*iter/gen_max);
        c2 =(0.5+2*iter/gen_max);
        bb = rand(ps,1);
        cc = rand(ps,1);
        bb = repmat(bb,1,D);
        cc = repmat(cc,1,D);


        vTmp1=(w.*vel + c1.*bb.*(pbest-pos) + c2.*cc.*(gbestrep-pos));
        posTmp1 = pos + vTmp1;

        vTmp2=(rand(ps,D).*vel+bb.*(pwin-pos)+c3*cc.*(center-pos));
        posTmp2 = pos + vTmp2;

        bnd = floor(ps*ratio);
        vel(1:bnd,:)=vTmp2(1:bnd,:);
        vel(bnd+1:ps,:)=vTmp1(bnd+1:ps,:);
        pos(1:bnd,:)=posTmp2(1:bnd,:);
        pos(bnd+1:ps,:)=posTmp1(bnd+1:ps,:);

        pos = ((pos>=Xmin)&(pos<=Xmax)).*pos...
            +(pos<Xmin).*(Xmin+0.25.*(Xmax-Xmin).*rand(ps,D))+...
            +(pos>Xmax).*(Xmax-0.25.*(Xmax-Xmin).*rand(ps,D));
        eval = feval(fhd,pos',fid);
        bin = (pbestval>eval)';
        if (sum(bin)>0)
            ratio = sum(bin(1:bnd))/(sum(bin(1:bnd))+sum(bin(bnd+1:ps)));
            ratio = max(0.1,ratio);
            ratio = min(0.9,ratio);
        end
        pbest(bin==1,:) =  pos(bin==1,:);
        pbestval(bin==1) = eval(bin==1);
        [gbestval,gbestid] = min(pbestval);
        nfe = nfe + ps;
        gbest = pbest(gbestid,:);
        nbest = gbest;
       %% Local Search
        var=zeros(ps,D);
        ave=mean(pos);
        for ii=1:ps
            %FIPS
            size = length(neighbor(ii,:));
            fi=(4.1/size).*rand(size,D);
            FIP=sum(fi.*pbest(neighbor(ii,:),:))./sum(fi);
            var(ii,:)=var(ii,:)+(pos(ii,:)-ave)+sum(fi).*(FIP-pos(ii,:));
        end
        sigma = sum(var)/ps;
        if rand<0.5
            competitor=nbest.*(cosd(90*(1-iter/gen_max)))+sigma;
        else
            competitor=nbest.*(sind(90*(1-iter/gen_max)))+sigma;
        end
        % Added by M. Omran to check boundaries

         competitor = ((competitor>=Xmin)&(competitor<=Xmax)).*competitor...
            +(competitor<Xmin).*(Xmin+0.25.*(Xmax-Xmin).*rand(1,D))+...
            +(competitor>Xmax).*(Xmax-0.25.*(Xmax-Xmin).*rand(1,D));
        evalSg = feval(fhd,competitor',fid);
        nfe=nfe+1;
        if evalSg<gbestval
            pbestval(gbestid) = evalSg;
            pbest(gbestid,:)=competitor;
            gbest=competitor;
        end
        gbestrep = repmat(gbest,ps,1);

        %% curves
        ccurve(old_nfe+1:nfe) = min(pbestval);
        % Determine the diversity of the archive
        prum = mean(pos);
        prum_mat=repmat(prum,ps,1);
        diam2 = (sum( sum( (pos-prum_mat).*(pos-prum_mat) ) ) / ps);
        diversity = sqrt(diam2);
        dcurve(old_nfe+1:nfe) = diversity;
        %%

        iter = iter+1;
    end

 %% curves
best_f = min(pbestval);
ccurve(nfe:nfe_max) = best_f;
dcurve(nfe:nfe_max) = diversity;
 %%

end
