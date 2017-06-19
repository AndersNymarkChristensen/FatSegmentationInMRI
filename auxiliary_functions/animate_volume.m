function animate_volume(V,p,direction,range,clim)
% CHANGES (IN RANGE) MADE AND NOT TESTED
% animate_volume(V,p,direction,range)
% V volume
% p pause, can be 'yes','no', 0 or >0
% direction, can be 'x','y','z','xy','xyz' etc, default 'xyz'
% range, can be a sequence, default full range

tic, fprintf('Animate volume ... ')

if nargin<2 || isempty(p)
    p = 'no';
end

if nargin<3 || isempty(direction)
    direction = 'xyz';
end

if nargin<4 || isempty(range)
    range_flag = false;
else
    range_flag = true;
end

if nargin<5
    % for equal color scaling
    clim = [min(V(:)), max(V(:))];
end

for d = 1:length(direction)
    cla, title(''), hold on, axis xy equal
    if direction(d)=='x'
        Vdir = permute(V,[2 3 1]);
        xlab = 'y'; ylab = 'z';
    elseif direction(d)=='y'
        Vdir = permute(V,[1 3 2]);
        xlab = 'x'; ylab = 'z';
    else % 'z'
        Vdir = V;
        xlab = 'x'; ylab = 'y';
    end
    
    if ~range_flag
        range = 1:size(Vdir,3);
    end
    
    for i = 1:numel(range)
        imagesc(Vdir(:,:,range(i))',clim) % transpose since I insist on having xyz volume
        title(['slice ',direction(d),' ',num2str(range(i)),'/',num2str(size(Vdir,3))])
        xlabel(xlab), ylabel(ylab)
        smart_pause(p)
    end
end

disp(['finished at ',num2str(toc)])

end



