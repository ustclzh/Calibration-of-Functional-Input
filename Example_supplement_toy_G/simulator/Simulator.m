function y=Simulator(f,data)
f=f*data.KL';
ql=-1;
ur=0;
n=length(f);
K=exp(f);
invk=1./K;
a=ql;
kinvmean=integrat_interval(invk,n);
b=ur-a*kinvmean;
y=integrat_interval(a.*invk)+b;
y=y(2:1:n-1);
end

function y=integrat_interval(k,s)
n=length(k);
T=[0,k(1:n-1)/2]+[0,k(2:n)/2];
y=cumsum(T)/(n-1);
if nargin==2
    y=y(s);
end
end

