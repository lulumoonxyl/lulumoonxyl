function [q,N]=quadrant_check(deg)
idx1=find(deg>=0&deg<90);
idx2=find(deg>=90&deg<=180);
idx3=find(deg>=(-180)&deg<(-90));
idx4=find(deg>=(-90)&deg<0);

q(idx1,1)=1;
q(idx2,1)=2;
q(idx3,1)=3;
q(idx4,1)=4;
[N,na,nb]=unique(q);

end