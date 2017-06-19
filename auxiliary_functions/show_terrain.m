function show_terrain(V,s,p,range,clim,direction)
%show_terrain_cut(V,S,p,range,clim,direction)
% V volume, a (X,Y,Z) matrix of intensities.
% s cut a (X,Y,K) matrix of heights, where K is number of cuts or a cell of
% such matrices
% p pause, can be 'yes','no', 0 or >0
% range of slices to be shown, can be a sequence, default full range
% clim, optional color scaling limits
% direction, can be 'x' or 'y',' default 'x'

if nargin<3 || isempty(p)
    p = 'no';
end

if nargin<5 || isempty(clim)
    % for equal color scaling
    clim = [min(V(:)), max(V(:))];
end

if nargin>5 && ~isempty(direction) && direction == 'x'
    V = permute(V,[3 2 1]);
    xlab = 'y';
else
    V = permute(V,[3 1 2]);
    direction = 'y'; % all input other than 'x' is 'y'
    xlab = 'x';
end

if nargin<4 || isempty(range)
    range = 1:size(V,3); 
end

if iscell(s)
    K = zeros(1,numel(s));
    for c = 1:numel(s)
        K(c) = size(s{c},3);
    end    
else
    K = size(s,3); % may also be 1;
end
col = jet(sum(K(:))); % a unique color for each cut

for i = range(:)'
    imagesc(V(:,:,i),clim), axis xy image, hold on  
    for c = 1:numel(K)
        if iscell(s)
            sc = s{c};
        else
            sc = s;
        end
        for k=1:K(c)
            sk = sc(:,:,k);
            if direction == 'x'
                sk = permute(sk,[2,1]);
            end
            plot(sk(:,i),'Color',col(sum(K(1:c-1))+k,:))        
        end
    end
    title(['slice ',direction,'=',num2str(i),'/',num2str(size(V,3))])
    xlabel(xlab), ylabel('z') 
    smart_pause(p)
end
end


