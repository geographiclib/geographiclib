n=8;
q=[1:n+1]';
x = sqrt(q);
xo = [0;x(1:n-1);0;-x(n-1:-1:1)];
xe = [x(1:n+1);x(n:-1:2)];
ind = [0:2*n-1]';
xso = [0;x(1:n);x(n-1:-1:1);0;-x(1:n-1);-x(n:-1:1)];
xse = [x(1:n);0;-x(n:-1:1);-x(2:n);0;x(n:-1:2)];
inds = [0:4*n-1]';
%figure(1);plot(ind,xo, ind, xe,inds/2,xso,inds/2,xse);

y=pi*[0:2*n-1]'/n;
yy=2*pi*[0:200]'/200;

co=-imag(fft(xo));
co=co(1:n)/n;                           % begin and end with 0
j=[0:n-1];
zo=sin(j.*y)*co;
%figure(2);plot(y,xo(1:2*n),'o-',y,zo,'x');

ce=real(fft(xe));
ce=ce(1:n+1)/n;
j=[0:n];
ce(1)=ce(1)/2;
ce(n+1) = ce(n+1)/2;
ze=cos(j.*y)*ce;
%figure(3);plot(y,xe(1:2*n),'o-',y,ze,'x');

ys=pi*[0:4*n-1]'/(2*n);
js=[0:2*n];
j1=[0:n-1];
cso=-imag(fft(xso));
cso=cso(2:2:end)/(2*n);
cso=cso(1:n);
zso=sin((2*j1+1).*yy)*cso;
%figure(4);plot(ys,xso(1:4*n),'o-',yy,zso,'x');

cse=real(fft(xse));
cse=cse(2:2:end)/(2*n);
cse=cse(1:n);
zse=cos((2*j1+1).*yy)*cse;
%figure(5);plot(ys,xse(1:4*n),'o-',yy,zse,'x');

xoa=xo(1:2:end);xob=xo(2:2:end);
xea=xe(1:2:end);xeb=xe(2:2:end);
xsoa=xso(1:2:end);xsob=xso(2:2:end);
xsea=xse(1:2:end);xseb=xse(2:2:end);

co=-imag(fft(xo));
co=co(1:n)/n;                           % begin and end with 0
coa=-imag(fft(xoa)); coa=coa(1:n/2)/(n/2);
cob=-imag(fft(xob).*exp(-i*[0:n-1]'*pi/n)); cob=cob(2:n/2+1)/(n/2);
coq=[0;coa(2:end)+cob(1:end-1);cob(end);-coa(end:-1:2)+cob(end-1:-1:1)]/2;

co=-imag(fft(xo));
co=co(2:n+1)/n;                           % begin and end with 0
coa=-imag(fft(xoa)); coa=coa(2:n/2+1)/(n/2);
cob=-imag(fft(xob).*exp(-i*[0:n-1]'*pi/n)); cob=cob(2:n/2+1)/(n/2);
coq=[0;coa(2:end)+cob(1:end-1);cob(end);-coa(end:-1:2)+cob(end-1:-1:1)]/2;
% coa ends with 0
coq=[coa+cob;-coa(end-1:-1:1)+cob(end-1:-1:1);0]/2;

ce=real(fft(xe));
ce=ce(1:n+1)/n;
cea=real(fft(xea));cea=cea(1:n/2+1)/(n/2);
ceb=real(fft(xeb).*exp(-i*[0:n-1]'*pi/n));ceb=ceb(1:n/2+1)/(n/2);
% ceb ends with 0
ceq=[cea+ceb;cea(end-1:-1:1)-ceb(end-1:-1:1)]/2;

cso=-imag(fft(xso));
cso=cso(2:2:end)/(2*n);
cso=cso(1:n);
csoa=-imag(fft(xsoa));csoa=csoa(2:2:end)/n;csoa=csoa(1:n/2);
csob=-imag(fft(xsob).*exp(-i*[0:2*n-1]'*pi/(2*n)));
csob=csob(2:2:end)/n;csob=csob(1:n/2);
csoq=[csoa+csob;-csoa(end:-1:1)+csob(end:-1:1)]/2;

cse=real(fft(xse));
cse=cse(2:2:end)/(2*n);
cse=cse(1:n);
csea=real(fft(xsea));csea=csea(2:2:end)/n;csea=csea(1:n/2);
cseb=real(fft(xseb).*exp(-i*[0:2*n-1]'*pi/(2*n)));
cseb=cseb(2:2:end)/n;cseb=cseb(1:n/2);
cseq=[csea+cseb;csea(end:-1:1)-cseb(end:-1:1)]/2;

return

xxob = [sqrt([1:n]),-sqrt([n:-1:1])]';
ccob =-imag(fft(xxob).*exp(-i*[0:2*n-1]'*pi/(2*n)));
ccob = ccob(1:n+1)/n;
return

n=16;
x = rand(2*n,1); z = fft(x);
za=fft(x(1:2:end)); zb=fft(x(2:2:end)).*exp(-i*[0:n-1]'*pi/n);
zq=[za+zb;za-zb];
d=z-zq;

y=z(1:n); y(1)=y(1)+i*z(n+1);
ya=za(1:n/2); ya(1)=ya(1)+i*za(n/2+1);
yb=zb(1:n/2); yb(1)=yb(1)+i*zb(n/2+1);
yq=[ya+yb;0;conj(ya(end:-1:2))-conj(yb(end:-1:2))];

f=exp(-i*[0:n-1]'*pi/n);
f1=[f(1:n/4);(1-i)/sqrt(2)];
f1=[f1;conj(f1(end-1:-1:1))*-i];
f1=[f1;-conj(f1(end-1:-1:2))];
