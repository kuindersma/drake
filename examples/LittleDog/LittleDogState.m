classdef LittleDogState < SingletonCoordinateFrame
  
  methods
    function obj=LittleDogState(r)
      typecheck(r,'TimeSteppingRigidBodyManipulator');
      manipStateFrame = r.getManipulator().getStateFrame();
      coordinates = manipStateFrame.coordinates;
      obj = obj@SingletonCoordinateFrame('LittleDogState',length(coordinates),'x',coordinates);
    end
  end
end
