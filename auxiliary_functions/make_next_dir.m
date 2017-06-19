function nextname = make_next_dir(dirname)
% makes a new directory dirname_x where x is a next free spot
% e.g. if map_1 and map_2 exist but not map_3, then it makes map_3
% vand@dtu.dk

i = 0;
success = 0;
message = 'nonempty';
while success~=1 || ~isempty(message)
    i = i+1;
    [success,message] = mkdir('.',[dirname,'_',num2str(i)]);
end
nextname = [dirname,'_',num2str(i)];