function [ghline,ghmarker]=plot_cov(cent,cov,style,alpha)
% Copied from plotcov_sam

N=40;
if(nargin<4), alpha=.05; end

if(length(cent)==1),
    cent=[cent,cent];
end

sc=sqrtm(cov);

a=2*pi*[0:N]'/N;
x=sqrt(chi2inv(1-alpha,2))*[cos(a),sin(a)]*sc';

x=[cent(1)+x(:,1),cent(2)+x(:,2)];
if(nargin<3)
    ghmarker=plot(cent(1),cent(2),'o');hold on
    ghline=plot(x(:,1),x(:,2));
else
    if isnumeric(style)
        ghmarker=plot(cent(1),cent(2),'o','color',style);hold on
        ghline=plot(x(:,1),x(:,2),'color',style);
    else
        ghmarker=plot(cent(1),cent(2),[style 'o']);hold on
        ghline=plot(x(:,1),x(:,2),style);
    end
end
