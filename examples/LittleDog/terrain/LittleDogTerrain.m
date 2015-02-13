classdef LittleDogTerrain < RigidBodyTerrain

  properties
    heightmaps;  % a collection of RigidBodyHeightMapTerrain objects
    x_positions;
  end
  
  methods
    function obj = LittleDogTerrain(terrain_files,num_ccw_rotations)
      % @param terrain_files a cell array of named terrain files.  
      %   the terrain will be layed out one after another in the x
      %   direction
      % @param num_ccw_rotations a vector with one element per terrain file containing
      %   the number of times to rotate each board 90 degrees in the counter-clockwise direction.  @default zeros
      
      obj = obj@RigidBodyTerrain();

      if ~iscell(terrain_files), terrain_files={terrain_files}; end
      if nargin<2, num_ccw_rotations = zeros(size(terrain_files)); end
      
      terrain_path = fileparts(mfilename('fullpath'));
      
      % life is simpler because the terrain files have all been
      % preprocessed
      N = length(terrain_files);
      
      x_position = 0;
      for i=1:N
        load(terrain_files{i});
        
        for j=1:num_ccw_rotations
          y=x;
          x=y;
          Z=fliplr(rot90(Z));
        end
        
        obj.heightmaps{i} = RigidBodyHeightMapTerrain(x+x_position,y,Z);

        obj.x_positions(i)=x_position;
        x_position = x_position+x(end);
      end
      obj.x_positions(end+1)=x_position;
    end
    
    function geom = getRigidBodyContactGeometry(obj)
      geom = {};
      for i=1:numel(obj.heightmaps)
        this_geom = obj.heightmaps{i}.getRigidBodyContactGeometry();
        if ~isempty(this_geom)
          geom = horzcat(geom,{this_geom});
        end
      end
    end

    function geom = getRigidBodyShapeGeometry(obj)
      geom = {};
      for i=1:numel(obj.heightmaps)
        this_geom = obj.heightmaps{i}.getRigidBodyShapeGeometry();
        if ~isempty(this_geom)
          geom = horzcat(geom,{this_geom});
        end
      end
    end
        
    function [z,normal] = getHeight(obj,xy)
      n = size(xy,1);
      z = 0*xy(1,:);
      normal = zeros(3,n);
      for i=1:length(obj.heightmaps)
        ind = xy(1,:)>=obj.x_positions(i) & xy(1,:)<obj.x_positions(i+1);
        [z(ind),normal(:,ind)]=getHeight(obj.heightmaps{i},xy(:,ind));
      end
    end
    
    function plotTerrain(obj)
      clf; hold on;
      for i=1:length(obj.heightmaps)
        plotTerrain(obj.heightmaps{i});
      end
    end
  end
  
  methods (Static=true)
    function terrainTest()
      terrain = LittleDogTerrain({'terrainA','terrainB','terrainFC'},[0,0,0]);
      figure(1);
      plotTerrain(terrain);
      
      figure(2);
      [x,y]=meshgrid(0:0.01:1.5,0:.01:.5);
      Z=getHeight(terrain,[x(:)';y(:)']);
      mesh(x,y,reshape(Z,size(x)));
      axis equal;
      
      options.floating = true;
      r = RigidBodyManipulator(fullfile('..','LittleDog.urdf'),options);
      r = setTerrain(r,terrain);
      r = compile(r);
      
      v = r.constructVisualizer();
      v.draw(0,[0;0;.4;zeros(100,1)]);
    end
  end
end
