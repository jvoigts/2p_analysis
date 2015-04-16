[relativeTime,absoluteTime]=extractTimesFromPrairieviewXML(directory);
%
% [relativeTime,absoluteTime]=extractTimesFromPrairieviewXML(directory);
%
% automatically locates xml name in directory - change this if more than
% one xml is in there
%

d=dir([directory,'/*.xml'])

% xmlread(fullfile(basedir,d(1).name)); % blows up

% instead, parse for instances of <Frame relativeTime="0"

try % just in case
    fclose(f);
end

f=fopen(fullfile(basedir,d(1).name),'r');
tline= ' ';
c=0;
absoluteTime=[];
relativeTime=[];
run=1;
while run
    tline = fgets(f);
    run=tline~=-1; % chekc if end is reached
    if run
        s=  sscanf(strtrim(tline), ['<Frame relativeTime="%f" absoluteTime="%f']);
        
        if numel(s)>0
            c=c+1;
            relativeTime(c)=s(1);
            absoluteTime(c)=s(2);
        end;
        if mod(c,1000)==0;
            fprintf('got %d timestamps\n',c);
        end;
    end;
end;
