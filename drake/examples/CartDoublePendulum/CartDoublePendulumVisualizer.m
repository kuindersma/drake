classdef CartDoublePendulumVisualizer < Visualizer
  
  properties
    l1=0.5;
    l2=0.5;
  end
  
  methods
    function obj = CartDoublePendulumVisualizer(plant)
      typecheck(plant,'CartDoublePendulumPlant');
      obj = obj@Visualizer(plant.getOutputFrame);
      obj.l1 = plant.l1;
      obj.l2 = plant.l2;
    end
    
    function draw(obj,t,x)
      l1=obj.l1;
      l2=obj.l2;
      persistent hFig base a1 a2 raarm1 raarm2 wb lwheel;
      if (isempty(hFig))
        hFig = sfigure(25);
        set(hFig,'DoubleBuffer', 'on');
        
        a1 = l1;
        a2 = l2;
        av = pi*[0:.05:1];
        theta = pi*[0:0.05:2];
        wb = .3; hb=.15;
        aw = .01;
        wheelr = 0.05;
        lwheel = [-wb/2 + wheelr*cos(theta); -hb-wheelr + wheelr*sin(theta)]';
        base = [wb*[1 -1 -1 1]; hb*[1 1 -1 -1]]';
        arm1 = [aw*cos(av-pi/2) -a1+aw*cos(av+pi/2)
          aw*sin(av-pi/2) aw*sin(av+pi/2)]';
        arm2 = [aw*cos(av-pi/2) -a2+aw*cos(av+pi/2)
          aw*sin(av-pi/2) aw*sin(av+pi/2)]';
        raarm1 = [(arm1(:,1).^2+arm1(:,2).^2).^.5, atan2(arm1(:,2),arm1(:,1))];
        raarm2 = [(arm2(:,1).^2+arm2(:,2).^2).^.5, atan2(arm2(:,2),arm2(:,1))];
      end
      
      sfigure(hFig); cla; hold on; view(0,90);
      patch(x(1)+base(:,1), base(:,2),0*base(:,1),'b','FaceColor',[.3 .6 .4])
      patch(x(1)+lwheel(:,1), lwheel(:,2), 0*lwheel(:,1),'k');
      patch(x(1)+wb+lwheel(:,1), lwheel(:,2), 0*lwheel(:,1),'k');
      patch(x(1)+raarm1(:,1).*sin(raarm1(:,2)+x(2)-pi),-raarm1(:,1).*cos(raarm1(:,2)+x(2)-pi), 1+0*raarm1(:,1),'r','FaceColor',[.9 .1 0]);
      patch(x(1)+a1*sin(x(2))+raarm2(:,1).*sin(raarm2(:,2)+x(3)-pi),-a1*cos(x(2))-raarm2(:,1).*cos(raarm2(:,2)+x(3)-pi), 1+0*raarm2(:,1),'r','FaceColor',[.9 .1 0]);
      plot3(x(1)+a1*sin(x(2)), -a1*cos(x(2)),1, 'ko',...
        'MarkerSize',14,'MarkerFaceColor','b');
      plot3(x(1)+a1*sin(x(2))+a2*sin(x(3)), -a1*cos(x(2))-a2*cos(x(3)),1, 'ko',...
        'MarkerSize',14,'MarkerFaceColor','b')
      plot3(x(1),0,1.5,'k.')
      set(gca,'XTick',[],'YTick',[])
      
      axis image;
      axis([-2.5 2.5 -2.5*l1 2.5*l1]);
    end
  end
end
