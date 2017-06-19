function surf_terrain(s)

hold on, axis vis3d off
view(3), camlight('right')
lighting phong, material shiny

if iscell(s)
    K = zeros(1,numel(s));
    for c = 1:numel(s)
        K(c) = size(s{c},3);
    end
else
    K = size(s,3); % may also be 1;
end
col = jet(sum(K(:))); % a unique color for each cut

alpha = 0.6;

for c = 1:numel(K)
    if iscell(s)
        sc = s{c};
    else
        sc = s;
    end    
    for k=1:K(c)
        sk = sc(:,:,k);
        surf(sk,'EdgeColor','none','FaceColor',col(sum(K(1:c-1))+k,:),...
            'FaceAlpha',alpha)
    end   
end

