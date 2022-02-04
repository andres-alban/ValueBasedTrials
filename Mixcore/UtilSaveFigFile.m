function [ rval ] =  UtilSaveFigFile( fignum, dirname, filename, modifier, extension )
% UtilSaveFigFile saves a file with the specifications of teh input
%


% if a directory is specified, and the directory does not exist, then create it
% also start to build up the file name for saving the file
if ~isempty(dirname)      
    tmp = dir(dirname);
    if ~length(tmp) %~isdir(dirname)
        mkdir(dirname);
    end
    fullpath = sprintf('%s/',dirname);
else
    fullpath = '';
end
fullpath = [fullpath filename];
if ~isempty(modifier)
    fullpath = [fullpath '_' modifier];
end
figurepath = [fullpath '.' extension];
%toggle = ['-d' extension]
saveas(fignum, figurepath);
print(fignum,figurepath,'-deps')

rval = true;

end
