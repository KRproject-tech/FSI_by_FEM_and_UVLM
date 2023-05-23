function h=ptristream(FlowP,y_pos,style)


N=length(FlowP);

h=zeros(1,N);

hold on;
for p=1:N    
    h(p)=plot3(FlowP(p).x,y_pos + 0*FlowP(p).x,FlowP(p).y,style);    
end
hold off

