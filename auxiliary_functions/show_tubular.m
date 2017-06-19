function show_tubular(V,S,p,range,clim)
%show_tubular_cut(V,S,p,range,clim)
% V volume, a (X,Y,Z) matrix of intensities.
% S cut (x and y coordinates), can be either a cell or a matrix:
%   - if cell, a K element cell of (N,2,Z) matrices
%   - if matrix, a (N,2,Z,K) matrix
%   where N is a number of cut points
%   2 is for x and y coordinate, and K is number of cuts.
% p pause, can be 'yes','no', 0 or >0
% range of slices to be shown, can be a sequence, default full range
% clim, optional color scaling limits

if nargin<3 || isempty(p)
    p = 'no';
end

if nargin<4 || isempty(range)
    range = 1:size(V,3);
end

if nargin<5
    % for equal color scaling
    clim = [min(V(:)), max(V(:))];
end

if iscell(S)
    K = zeros(1,numel(S));
    for c = 1:numel(S)
        K(c) = size(S{c},4);
    end
else
    K = size(S,4); % may also be 1;
end
col = jet(sum(K(:))); % a unique color for each cut

for i = range(:)'
    imagesc(V(:,:,i),clim), colormap gray, hold on
    for c = 1:numel(K)
        if iscell(S)
            Sc = S{c};
        else
            Sc = S;
        end
        for k = 1:K(c)
            Sk = Sc(:,:,:,k);
            plot(Sk([1:end,1],i,1),Sk([1:end,1],i,2),'Color',col(sum(K(1:c-1))+k,:))
        end
    end    
    title(['slice z ',num2str(i),'/',num2str(size(V,3))])
    xlabel('y'), ylabel('x'), axis xy image
%     saveas(gcf, ['..\data_examples\goactiwe\Figures\tubular',num2str(i),'.png'])
    smart_pause(p)
end

end



