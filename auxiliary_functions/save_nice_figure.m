function save_nice_figure(name,lw,fs,leave_handles,do_save)
% sets line width and font size to 3 and 15 (or other values) but will
% leave objects from leave_handles list unchanged
% saves current figure under name (if name not empty and do_save not false)
% vand@dtu.dk

if nargin<5
    do_save = true;
end
if nargin<4
    leave_handles=[];
end
if nargin<3 || isempty(fs)
    fs = 15;
end
if nargin<2 || isempty(lw)
    lw = 3;
end
if isempty(lw), lw=3; end
if isempty(fs), fs=15; end

h_obj = findobj(gcf); % all non-hidden handles in current figure

% set(gca,'LineWidth',lw,'FontSize',fs)
% findobj(gcf,'Type','axes');
% set(findall(gca,'Type','line'),'LineWidth',lw)
% set(findobj(gcf,'Type','line'),'LineWidth',lw); %
% set(findall(gca,'Type','text'),'FontSize',fs)

set(setdiff(findobj(h_obj,'Type','line'),leave_handles),'LineWidth',lw)
% findall finds also hidden handles, so we get to title and labels
set(setdiff(findall(h_obj,'Type','text'),leave_handles),'FontSize',fs)
set(setdiff(findobj(h_obj,'Type','axes'),leave_handles),'LineWidth',lw,'FontSize',fs)

box on

if isempty(name)
    do_save = false;
elseif numel(name)<4 || name(end-3)~='.'
    name = [name,'.png']; 
end
if do_save
    saveas(gcf,name)
end