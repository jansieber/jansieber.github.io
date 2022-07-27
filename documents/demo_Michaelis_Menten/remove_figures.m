function namelist=remove_figures(varstruct)
%% Save all variables that are not graphics
% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics
% object 
%%
% from
% <https://stackoverflow.com/questions/38131166/save-matlab-workspace-without-saving-or-deleting-figures>
tosave = cellfun(@isempty,regexp({varstruct.class},'^matlab\.(ui|graphics)\.'));
namelist={varstruct(tosave).name};
end
