function [xlist,ylist]=adaptive(func,opt)
if ~isa(func,'function_handle')
    error('%s is not a function handle',func2str(func));
end
if ~isa(opt,'struct')
    error('option is not a struct');
end
if isfield(opt,{'min','max'})
    if opt.min>=opt.max
        error('min of plot range %f is larger than max %f',opt.min,opt.max);
    end
else
    error('plot range is missing');
end
if ~isfield(opt,'points')
    opt.points=50;
end
if ~isfield(opt,'maxrecursion')
    opt.maxrecursion=5;
end
if ~isfield(opt,'angle')
    opt.angle=5/360*2*pi;
end
xlist=linspace(opt.min,opt.max,opt.points);

dx=(opt.max-opt.min)/(opt.points-1);
ylist=zeros(1,opt.points);
parfor i=1:opt.points
    ylist(i)=func(xlist(i));
end
dylist=diff(ylist);
l1sq=dx^2+dylist(1:end-1).^2;
l2sq=dx^2+dylist(2:end).^2;
queue=(l1sq+l2sq-2*sqrt(l1sq.*l2sq)*cos(opt.angle)>=diff(dylist).^2);
% xmin=opt.min;
% xmax=opt.max;
for i=1:length(queue)
    if queue(i)==0
        % for each i-th point, if the second deriviative exceeds threshold
        % angle, we just refine it to grid on [i-1/2,i+1/2]. if they are on
        % the boundary, [min,i+1/2] or [i-1/2,max] will be adopted.
        opt2=opt;
        if i~=1
            opt2.min=(opt.max-opt.min)/(length(queue)+1)*(i-1/2)+opt.min;
        end
        if i~=length(queue)
            opt2.max=(opt.max-opt.min)/(length(queue)+1)*(i+1/2)+opt.min;
        end
        [xlist2,ylist2]=adaptive(func,opt2);
        xlist=[xlist,xlist2];
        ylist=[ylist,ylist2];
    end
end 
end
    