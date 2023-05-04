function q=check_qua_turn(deg)

idx1=find(deg>=(-45)&deg<45);
idx2=find(deg>=45&deg<135);
idx3=find((deg>=135&deg<=180)|(deg>=(-180)&deg<(-135)));
idx4=find(deg>=(-135)&deg<(-45));

q(idx1,1)=1;
q(idx2,1)=2;
q(idx3,1)=3;
q(idx4,1)=4;


end 