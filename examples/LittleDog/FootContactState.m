classdef FootContactState < CoordinateFrame & Singleton
  methods
    function obj=FootContactState()
      obj = obj@CoordinateFrame('FootContactState',4,'x',{'right_front','left_front','right_back','left_back',});
      obj = obj@Singleton();
    end
  end
end
