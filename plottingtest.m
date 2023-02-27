figure;

plt1 = subplot(4,1,1);
plot(sensors.time(60000:70000),F_feet(60000:70000,3,1));
title('RF','FontSize',16);
plt2 = subplot(4,1,2);
plot(sensors.time(60000:70000),F_feet(60000:70000,3,2));
title('LF','FontSize',16);
plt3 = subplot(4,1,3);
plot(sensors.time(60000:70000),F_feet(60000:70000,3,3));
title('RH','FontSize',16);
plt4 = subplot(4,1,4);
plot(sensors.time(60000:70000),F_feet(60000:70000,3,4));
title('LH','FontSize',16);
xlabel('Time (s)','FontSize',16);

p1 = get(plt1,'position');
p2 = get(plt2,'position');
p3 = get(plt3,'position');
p4 = get(plt4,'position');
height = p1(2) + p1(4) + p2(4) - p2(2) - p3(2);
h3 = axes('position',[(p3(1)+p2(1))/2-0.03 (p3(2)+p2(2))/2 p2(3) height],'visible','off');
h_label = ylabel('Ground Reaction Force (N)','visible','on','FontSize',16);
