
close all;clear all;clc;
load rotUni.mat
f = figure(1);
hold on;
filename = 'rotUni.gif';
xlabel('x')
ylabel('y')
zlabel('z')
axis([-1.1 1.1 -1.1 1.1 -.3 1.1])
view([-135 30])
grid on
v = VideoWriter('rotUni.avi');
open(v);

qsT = qs;
ii = 1;
while ii < size(qsT,2)
    qsT(:,ii) = [];
    ii = ii + 2;
end

initVec = q2Q(qs(:,1))*eye(3);
finalVec = q2Q(qs(:,end))*eye(3);
for ii = 1:size(qs,2)
    if ii <= size(qsT,2)
        Q1 = q2Q(qsT(:,ii));
        vecs1 = Q1*eye(3);
        quiver3(0,0,0,vecs1(1,1),vecs1(1,2),vecs1(1,3),'r','LineWidth',1);
        quiver3(0,0,0,vecs1(2,1),vecs1(2,2),vecs1(2,3),'r','LineWidth',1);
        quiver3(0,0,0,vecs1(3,1),vecs1(3,2),vecs1(3,3),'r','LineWidth',4);
    else
        quiver3(0,0,0,finalVec(1,1),finalVec(1,2),finalVec(1,3),'r','LineWidth',1);
        quiver3(0,0,0,finalVec(2,1),finalVec(2,2),finalVec(2,3),'r','LineWidth',1);
        quiver3(0,0,0,finalVec(3,1),finalVec(3,2),finalVec(3,3),'r','LineWidth',4);
    end
    Q2 = q2Q(qOpt(:,ii));
    vecs2 = Q2*eye(3);
    quiver3(0,0,0,vecs2(1,1),vecs2(1,2),vecs2(1,3),'b','LineWidth',1);
    quiver3(0,0,0,vecs2(2,1),vecs2(2,2),vecs2(2,3),'b','LineWidth',1);
    quiver3(0,0,0,vecs2(3,1),vecs2(3,2),vecs2(3,3),'b','LineWidth',4);
    pause(dts(1)/1000);
    if ii ~= size(qs,2);
%         frame = getframe(f);
%         im = frame2im(frame);
%         [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
%         if ii == 1
% %             imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%             writeVideo(v,frame)
%         else
%             writeVideo(v,frame)
%         end
        cla
    end
end
close(v);

% figure;
% axis([-1.1 1.1 -1.1 1.1 -.3 1.1])
% hold on
%
% initVec = quat2Q(qOpt(:,1))*eye(3);
% for ii = 1:size(qs,2)
%     Q = quat2Q(qOpt(:,ii));
%     vecs = Q*eye(3);
%     quiver3(zeros(3,1), zeros(3,1), zeros(3,1), vecs2(:,1), vecs2(:,2), vecs2(:,3),'b','LineWidth',2);
%     pause(dt/10);
%     if ii ~= size(qs,2);
%         cla
%         quiver3(zeros(3,1), zeros(3,1), zeros(3,1), initVec(:,1), initVec(:,2), initVec(:,3),'g','LineWidth',2);
%     end
% end