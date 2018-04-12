function no = BARON_uninstall(del,removeAll)
%BARON_UNINSTALL  BARON Uninstaller
%
% Simply run this function to remove all instances of the MATLAB-BARON
% (matbar) interface on your computer. Optionally, call with del=1 to also
% delete all matbar interfaces.


%Token we search for to find BARON interface directories
token = 'BARON_install.m';
no = 0;
%Check nargin in, default don't delete
if(nargin < 1 || isempty(del))
    del = 0;
end
if(nargin < 2 || isempty(removeAll))
    removeAll = 1; %for removing all versions (including current)
end
%Check if we have anything to remove
paths = which(token,'-all');
len = length(paths);
if(~len) %should always be at least 1 if we are in correct directory
    error('Expected to find "%s" in the current directory - please ensure you are in the BARON Interface (matbar) directory');        
elseif(len == 1 && removeAll==0)
    %if len == 1, either we are in the correct folder with nothing to remove, or we are in the
    %wrong folder and there are files to remove, check CD
    if(any(strfind(paths{1},cd)))
        no = 0;
        return;
    else
        error('Expected to find "%s" in the current directory - please ensure you are in the BARON Interface (matbar) directory');
    end    
else %old ones to remove
    %Remove each folder found, and all subdirs under
    for n = 2-removeAll:len
        %Absolute path to remove
        removeP = paths{n};
        %Search backwards for first file separator (we don't want the filename)
        for j = length(removeP):-1:1
            if(removeP(j) == filesep)
                break;
            end
        end
        removeP = removeP(1:max(j-1,1));        

        %Everything is lowercase to aid matching
        lrpath = lower(removeP);
        opath = regexp(lower(path),';','split');

        %Find & Remove Matching Paths
        no = 0;
        for i = 1:length(opath)
            %If we find it in the current path string, remove it
            fnd = strfind(opath{i},lrpath);        
            if(~isempty(fnd))  
                rmpath(opath{i});
                no = no + 1;
            end
        end
        rehash;
        %If delete is specified, also delete the directory
        if(del)
            stat = recycle; recycle('on'); %turn on recycling
            try
                rmdir(removeP,'s'); %not sure if we dont have permissions here
            catch
                %just do the best we can, we can't delete the folder this
                %function is running in!
            end
            recycle(stat); %restore to original
        end
    end 
    %If we found something, save path changes
    if(no > 0)
        try
            savepath;
        catch %#ok<CTCH>
            warning('baron:uninstall',['It appears you do not have administrator rights on your computer to save the Matlab path. '...
                                       'In order to uninstall BARON you will need to contact your system administrator to obtain administrator rights.']);
            return; %don't delete - too many warnings
        end
        if(nargout < 1)
            fprintf('MATLAB-BARON Interface successfully uninstalled.\n');
        end
    end
end