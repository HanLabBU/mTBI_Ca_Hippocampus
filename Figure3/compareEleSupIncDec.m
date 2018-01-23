%Takes the output from activityChanges.m and plotMeanValuesPopulation.m and
%compares them

%expmeans/ctlmeans is mean elevated/decreasing values of population
%exppcts/ctlpcts is for increasing/decreasing activity cells

activity = '\\engnas.ad.bu.edu\Research\eng_research_handata\Kyle\Code\TBI_Paper_Scripts\Figure3\activityResults_20170918_Day2.mat';
baseline = '\\engnas.ad.bu.edu\Research\eng_research_handata\Kyle\Code\TBI_Paper_Scripts\Figure3\sigpopulation_20170918_Day2.mat';

load(activity)
load(baseline)

experimental = struct();
control = struct();
for idx = 1:numel(expmeans)
    %Experimental
    einc = exppcts(idx).inc;
    edec = exppcts(idx).dec;
    enm = ~(einc+edec);
    eup = expmeans(idx).up;
    edown = expmeans(idx).down;
    econ = ~(eup+edown);
    
    %Experimental
    eupinc = zeros(numel(eup),1);
    eupdec = zeros(numel(eup),1);
    eupnone = zeros(numel(eup),1);
    edowninc = zeros(numel(edown),1);
    edowndec = zeros(numel(edown),1);
    edownnone = zeros(numel(edown),1);
    econinc = zeros(numel(econ),1);
    econdec = zeros(numel(econ),1);
    econnone = zeros(numel(econ),1);
    %Comparisons
    for iter = 1:numel(eupinc)
        if eup(iter) == 1
            if einc(iter) == 1
                eupinc(iter) = 1;
            elseif edec(iter) == 1
                eupdec(iter) = 1;
            elseif enm(iter) == 1
                eupnone(iter) = 1;
            end
        end
        if edown(iter) == 1
            if einc(iter) == 1
                edowninc(iter) = 1;
            elseif edec(iter) == 1
                edowndec(iter) = 1;
            elseif enm(iter) == 1
                edownnone(iter) = 1;
            end
        end
        if econ(iter) == 1
            if einc(iter) == 1
                econinc(iter) = 1;
            elseif edec(iter) == 1
                econdec(iter) = 1;
            elseif enm(iter) == 1
                econnone(iter) = 1;
            end
        end
    end
    experimental(idx).upinc = eupinc;
    experimental(idx).updec = eupdec;
    experimental(idx).upnone = eupnone;
    experimental(idx).downinc = edowninc;
    experimental(idx).downdec = edowndec;
    experimental(idx).downnone = edownnone;
    experimental(idx).coninc = econinc;
    experimental(idx).condec = econdec;
    experimental(idx).connone = econnone;
    experimental(idx).tot = numel(eupinc);
end
for idx=1:numel(ctlmeans)
    %Control
    cinc = ctlpcts(idx).inc;
    cdec = ctlpcts(idx).dec;
    cnm = ~(cinc+cdec);
    cup = ctlmeans(idx).up;
    cdown = ctlmeans(idx).down;
    ccon = ~(cup+cdown);
    
    %Control
    cupinc = zeros(numel(cup),1);
    cupdec = zeros(numel(cup),1);
    cupnone = zeros(numel(cup),1);
    cdowninc = zeros(numel(cdown),1);
    cdowndec = zeros(numel(cdown),1);
    cdownnone = zeros(numel(cdown),1);
    cconinc = zeros(numel(ccon),1);
    ccondec = zeros(numel(ccon),1);
    cconnone = zeros(numel(ccon),1);
    %Comparisons
    for iter = 1:numel(cupinc)
        if cup(iter) == 1
            if cinc(iter) == 1
                cupinc(iter) = 1;
            elseif cdec(iter) == 1
                cupdec(iter) = 1;
            elseif cnm(iter) == 1
                cupnone(iter) = 1;
            end
        end
        if cdown(iter) == 1
            if cinc(iter) == 1
                cdowninc(iter) = 1;
            elseif cdec(iter) == 1
                cdowndec(iter) = 1;
            elseif cnm(iter) == 1
                cdownnone(iter) = 1;
            end
        end
        if ccon(iter) == 1
            if cinc(iter) == 1
                cconinc(iter) = 1;
            elseif cdec(iter) == 1
                ccondec(iter) = 1;
            elseif cnm(iter) == 1
                cconnone(iter) = 1;
            end
        end
    end
    control(idx).upinc = cupinc;
    control(idx).updec = cupdec;
    control(idx).upnone = cupnone;
    control(idx).downinc = cdowninc;
    control(idx).downdec = cdowndec;
    control(idx).downnone = cdownnone;
    control(idx).coninc = cconinc;
    control(idx).condec = ccondec;
    control(idx).connone = cconnone;
    control(idx).tot = numel(cupinc);
end

%Totals
expall = struct();
ctlall = struct();
for idx = 1:numel(experimental);
    %Experimental
    expall.upinc(idx) = sum(experimental(idx).upinc);
    expall.updec(idx) = sum(experimental(idx).updec);
    expall.upnone(idx) = sum(experimental(idx).upnone);
    expall.downinc(idx) = sum(experimental(idx).downinc);
    expall.downdec(idx) = sum(experimental(idx).downdec);
    expall.downnone(idx) = sum(experimental(idx).downnone);
    expall.coninc(idx) = sum(experimental(idx).coninc);
    expall.condec(idx) = sum(experimental(idx).condec);
    expall.connone(idx) = sum(experimental(idx).connone);
    expall.upnum(idx) = sum([expall.upinc(idx),expall.updec(idx),expall.upnone(idx)]);
    expall.downnum(idx) = sum([expall.downinc(idx),expall.downdec(idx),expall.downnone(idx)]);
    expall.connum(idx) = sum([expall.coninc(idx),expall.condec(idx),expall.connone(idx)]);
    %Percentages
    expall.upincpct(idx) = sum(experimental(idx).upinc)/expall.upnum(idx);
    expall.updecpct(idx) = sum(experimental(idx).updec)/expall.upnum(idx);
    expall.upnonepct(idx) = sum(experimental(idx).upnone)/expall.upnum(idx);
    expall.downincpct(idx) = sum(experimental(idx).downinc)/expall.downnum(idx);
    expall.downdecpct(idx) = sum(experimental(idx).downdec)/expall.downnum(idx);
    expall.downnonepct(idx) = sum(experimental(idx).downnone)/expall.downnum(idx);
    expall.conincpct(idx) = sum(experimental(idx).coninc)/expall.connum(idx);
    expall.condecpct(idx) = sum(experimental(idx).condec)/expall.connum(idx);
    expall.connonepct(idx) = sum(experimental(idx).connone)/expall.connum(idx);
end
for idx = 1:numel(control)
    %Control
    ctlall.upinc(idx) = sum(control(idx).upinc);
    ctlall.updec(idx) = sum(control(idx).updec);
    ctlall.upnone(idx) = sum(control(idx).upnone);
    ctlall.downinc(idx) = sum(control(idx).downinc);
    ctlall.downdec(idx) = sum(control(idx).downdec);
    ctlall.downnone(idx) = sum(control(idx).downnone);
    ctlall.coninc(idx) = sum(control(idx).coninc);
    ctlall.condec(idx) = sum(control(idx).condec);
    ctlall.connone(idx) = sum(control(idx).connone);
    ctlall.upnum(idx) = sum([ctlall.upinc(idx),ctlall.updec(idx),ctlall.upnone(idx)]);
    ctlall.downnum(idx) = sum([ctlall.downinc(idx),ctlall.downdec(idx),ctlall.downnone(idx)]);
    ctlall.connum(idx) = sum([ctlall.coninc(idx),ctlall.condec(idx),ctlall.connone(idx)]);
    %Percentages
    ctlall.upincpct(idx) = sum(control(idx).upinc)/ctlall.upnum(idx);
    ctlall.updecpct(idx) = sum(control(idx).updec)/ctlall.upnum(idx);
    ctlall.upnonepct(idx) = sum(control(idx).upnone)/ctlall.upnum(idx);
    ctlall.downincpct(idx) = sum(control(idx).downinc)/ctlall.downnum(idx);
    ctlall.downdecpct(idx) = sum(control(idx).downdec)/ctlall.downnum(idx);
    ctlall.downnonepct(idx) = sum(control(idx).downnone)/ctlall.downnum(idx);
    ctlall.conincpct(idx) = sum(control(idx).coninc)/ctlall.connum(idx);
    ctlall.condecpct(idx) = sum(control(idx).condec)/ctlall.connum(idx);
    ctlall.connonepct(idx) = sum(control(idx).connone)/ctlall.connum(idx);
end

[pupinc, hupinc, statsupinc] = ranksum(ctlall.upincpct, expall.upincpct);
[pupdec, hupdec, statsupdec] = ranksum(ctlall.updecpct, expall.updecpct);
[pupnone, hupnone, statsupnone] = ranksum(ctlall.upnonepct, expall.upnonepct);
[pdowninc, hdowninc, statsdowninc] = ranksum(ctlall.downincpct, expall.downincpct);
[pdowndec, hdowndec, statsdowndec] = ranksum(ctlall.downdecpct, expall.downdecpct);
[pdownnone, hdownnone, statsdownnone] = ranksum(ctlall.downnonepct, expall.downnonepct);