function [rotmat,names]=get_rotation(varargin)
default={'rotmat',eye(2),'names',{'x','y'}};
options=ScSetOptions(default,varargin,'pass_on');
rotmat=options.rotmat;
names=options.names;
end