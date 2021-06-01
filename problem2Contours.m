function problem2Contours(edges)
    if nargin<1
        edges = 5;
    end

    Himmelblau = @(x1,x2)(sin(x1).^2+x2-11).*cos(x2.*x2) + (x1+cos(x2).^2-7).*sin(x1);

    % Compute the function values
    x1 = linspace(-edges,edges, 400); % linspace(-5,5, 100);
    x2 = linspace(-edges,edges, 400); % linspace(-5,5, 100);
    [X1,X2]=meshgrid(x1,x2);
    F = Himmelblau(X1,X2);
    

    % Make contour plot
    if edges == 5
        v = [-20:5:20];
    else
        v = [-20:5:20 100:20:200 200:floor((max(max(F))-200)/(max(max(F))/200)):max(F) max(max(F))];
    end
    hold on
    [C,h] = contour(X1,X2,F,v,'linewidth',2);
    


    set( get( get( h, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    xlabel('x_1','Fontsize',14)
    ylabel('x_2','Fontsize',14)
    set(gca,'fontsize',14);
    set(gcf, 'Color', 'w');
    colorbar
    axis image