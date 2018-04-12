function BARON_install
% Installation File for BARON

% In order to run this tool, please run this file to setup the required
% directories. You MUST be in the current directory of this file!

cpath = cd;
ofile = which('pathdef.m');

try
    cd('Interface');
catch %#ok<CTCH>
    error('You don''t appear to be in the BARON Interface (matbar) directory');
end
cur_ver = baron();
cd(cpath);

fprintf('\n------------------------------------------------\n')
fprintf(['  INSTALLING BARON INTERFACE ' cur_ver '\n'])

fprintf('\n- Checking for previous versions of the MATLAB-BARON (matbar) Interface...\n');
%token file to search for, must be in the base directory of the toolbox
no = BARON_uninstall(0,1); %change 0 to 1 to delete old versions
if(no < 1)
    fprintf('Could not find a previous installation.\n');
else
    fprintf('Successfully uninstalled previous version(s).\n');
end

fprintf('\n- Adding BARON Paths to MATLAB Search Path...');
genp = genpath(cd);
genp = regexp(genp,';','split');
%Folders to exclude from adding to Matlab path
i = 1;
rInd{:,:,i} = strfind(genp,'vti_cnf'); i = i + 1;
rInd{:,:,i} = strfind(genp,'vti_pvt'); i = i + 1;
rInd{:,:,i} = strfind(genp,'.git'); i = i + 1;
rInd{:,:,i} = strfind(genp,'tex'); i = i + 1;
ind = NaN(length(rInd{1}),1);
%Track indices of paths to remove from list
for i = 1:length(rInd{1})
    for j = 1:size(rInd,3)
        if(any(rInd{j}{i}))
            ind(i) = 1;
        end
    end
end
%Remove paths from above and add to matlab path
genp(ind == 1) = [];
addpath(genp{:});
rehash


%Check we can write to pathdef (linux issue)
fid = fopen(ofile,'w'); istemp = false; isok = true;
if(fid < 0)
    fprintf('\n\n'); istemp = true;
    switch(computer)
        case {'PCWIN','PCWIN64'}
            warning('baron:install',['It appears you do not have administrator rights on your computer to save the Matlab path. '...
                             'In order to run BARON you will need to install it each time you wish to use it. To fix '...
                             'this please contact your system administrator to obtain administrator rights.']);
        otherwise
            warning('baron:install',['Please run MATLAB as Super User (sudo) in order to install the BARON Interface paths.\n'...
                                     'Currently you will need to install BARON each time you wish to use it!']);
    end
else
    try
        fclose(fid);
        savepath;
        fprintf('Done\n');
    catch %#ok<CTCH>
        warning('baron:install',['It appears you do not have administrator rights on your computer to save the Matlab path. '...
                                 'In order to run BARON you will need to install it each time you wish to use it. To fix '...
                                 'this please contact your system administrator to obtain administrator rights.']);
    end
end

%Check for executable
fprintf('\n- Checking for the BARON executable (barin)...');
barloc = baron('exe');
if(isempty(barloc))
    isok = false;
    intloc = which('baron.m'); idx = find(intloc==filesep); intloc = intloc(1:idx(end)); 
    err = sprintf(['\nCannot find the BARON executable! Please complete the following steps:\n1) Download the executable for your operating system from: '...
                   'http://minlp.com/download\n2) Rename the executable "barin"\n3) Place it in the following folder: '...
                   '%s\n4) Re-run this installer.\n'],intloc);
    switch(computer)
        case{'PCWIN','PCWIN64'}
            err = regexprep(err,filesep,[filesep filesep filesep filesep]);
    end
    fprintf(2,err);
else
    fprintf('Done\n');
end

if(istemp)
    fprintf('\nBARON Interface Temporarily Installed.\n');
elseif(isok)
    fprintf('\n  BARON Interface Installation Complete!\n');
else
    fprintf(2,'\n  BARON Interface Installation Failed.\n');
end
    

%Clean up
clear cpath cur_ver genp fid ofile ans

disp('------------------------------------------------')