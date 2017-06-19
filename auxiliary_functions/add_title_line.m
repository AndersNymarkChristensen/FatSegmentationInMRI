function add_title_line(title_line,axis_handle)
% adds a line on top of the current title
% vand@dtu.dk, 17.03.2014

if nargin<2
    axis_handle = gca;
end

title_text =get(get(axis_handle,'Title'),'String');
if iscell(title_text)
    title_text(2:end+1) = title_text(1:end);
    title_text{1} = title_line;
else
    title_text = {title_line,title_text};
end
title(title_text)