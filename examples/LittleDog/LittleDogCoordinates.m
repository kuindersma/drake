classdef LittleDogCoordinates < SingletonCoordinateFrame
  % atlas q
  methods
    function obj=LittleDogCoordinates(r)
      typecheck(r,'TimeSteppingRigidBodyManipulator');
      nq = r.getNumPositions();
      obj = obj@SingletonCoordinateFrame('LittleDogCoordinates',nq,'x',r.getStateFrame.coordinates(1:nq));
    end
  end
end
