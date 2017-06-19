function smart_pause(p)
% p can be 'yes','no', 0 or >0
% vand@dtu.dk

if strcmp(p,'yes')
    pause
elseif strcmp(p,'no')||p==0
    drawnow
else
    pause(p)
end

end