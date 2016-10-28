function [maxval,new_atom]=lmo_ellipse(s,param)
if size(s,2)>1,
    error('s is not a column vector');
end
if size(s,1)~=2,
    error('This code is for dimension 2 only is not a column vector');
end

a=param.a;
b=param.b;

mu=sqrt((s(1)^2)/a+(s(2)^2)/b);
new_atom=[s(1)/a;s(2)/b]./mu;

maxval=new_atom'*s;

if 0
    p0=s(2)*b/(s(1)*a);
    p=mod(p0,pi)-pi/2;
    new_atom01=[a*cos(p);b*sin(p)];
    maxval01=new_atom01'*s;
    new_atom02=[a*cos(p+pi/2);b*sin(p+pi/2)];
    maxval02=new_atom02'*s;
    new_atom03=[a*cos(p+pi);b*sin(p+pi)];
    maxval03=new_atom03'*s;
    new_atom04=[a*cos(p+pi*3/2);b*sin(p+pi*3/2)];
    maxval04=new_atom04'*s; 
    
    t=linspace(-pi,pi,100);
    x=a*cos(t);
    y=b*sin(t);
    figure(1);clf;
    plot(x,y, '.'); hold on;
    plot([0 s(1)], [0 s(2)],'b-','lineWidth',2); hold on;
    plot([0 new_atom(1)], [0 new_atom(2)],'r-','lineWidth',2); hold on;
    plot([0 new_atom01(1)], [0 new_atom01(2)],'g-','lineWidth',2); hold on;
    plot([0 new_atom02(1)], [0 new_atom02(2)],'m-','lineWidth',2); hold on;
    plot([0 new_atom03(1)], [0 new_atom03(2)],'k-','lineWidth',2); hold on;
    plot([0 new_atom04(1)], [0 new_atom04(2)],'y-','lineWidth',2); hold on;
    axis equal
    keyboard;
end

end
