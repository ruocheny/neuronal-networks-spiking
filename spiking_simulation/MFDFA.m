% code source: Ihlen, Espen AF. "Introduction to multifractal detrended fluctuation 
% analysis in Matlab." Frontiers in physiology 3 (2012): 141.

function [Fq, Hq, qRegLine, scale, tq, hq, Dq] = MFDFA(ts)

q=-5:0.25:5;
X=cumsum(ts-mean(ts));

if length(ts) <= 64 % skip shorter time series
    Fq = nan;
    Hq = nan;
    qRegLine = nan;
    scale = nan;
    tq = nan;
    hq = nan;
    Dq = nan;
    
else
    smin = 8;
    smax = floor(length(ts)/4);
    scale = log2(smin):1:log2(smax);
    scale = 2.^scale;
    m=1;
    for ns=1:length(scale)
        segments(ns)=floor(length(X)/scale(ns));
        for v=1:segments(ns)
            Index=((((v-1)*scale(ns))+1):(v*scale(ns)));
            C=polyfit(Index,X(Index),m);
            fit=polyval(C,Index);
            RMS{ns}(v)=sqrt(mean((X(Index)-fit).^2));
        end
        for nq=1:length(q)
            qRMS{nq,ns}=RMS{ns}.^q(nq);
            Fq(nq,ns)=mean(qRMS{nq,ns}).^(1/q(nq));
        end
        Fq(q==0,ns)=exp(0.5*mean(log(RMS{ns}.^2)));
    end
    
    for nq=1:length(q)
        C=polyfit(log2(scale),log2(Fq(nq,:)),1);
        Hq(nq)=C(1);
        qRegLine{nq}=polyval(C,log2(scale));
    end
    
    tq=Hq.*q-1; % q-order mass exponent (tau)
    hq=diff(tq)./(q(2)-q(1)); % alpha
    Dq=(q(1:end-1).*hq)-tq(1:end-1); %f(alpha)
    
end