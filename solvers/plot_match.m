function plot_match(a, X, inliers, plotOutliers, numPlot)

    N = size(X, 2); 
    outliers = 1:N; 
    outliers(inliers) = []; 

    im1 = a.im1;
    im2 = a.im2;
    
    [rows, columns, numberOfColorChannels] = size(im1);
    [rows2, columns2, numberOfColorChannels2] = size(im2);
    if(numberOfColorChannels == 3 && numberOfColorChannels2 == 1)
        im1 = im2single(rgb2gray(im1));
    elseif (numberOfColorChannels == 1 && numberOfColorChannels2 == 3)
        im2 = im2single(rgb2gray(im2));
    end


    dh1 = max(size(im2,1)-size(im1,1),0);
    dh2 = max(size(im1,1)-size(im2,1),0);

    if (numPlot > 0)
        inliers = inliers(1:min(numPlot,length(inliers)));
        outliers = outliers(1:min(numPlot, length(outliers)));
    end
    

    %subplot(1,1,1)
    subtightplot(1,1,1,[0.01,0.01],0.01,0.01)
    imagesc([padarray(im1,dh1,'post') padarray(im2,dh2,'post')]); hold on;
    o = size(im1,2);
    if (plotOutliers<2)
    plot(a.X1(1, inliers),a.X1(2, inliers),'+b', 'MarkerSize',4, 'linewidth', 1.2); 
    plot(a.X2(1, inliers)+o,a.X2(2, inliers),'+b', 'MarkerSize',4, 'linewidth', 1.2);    
    line([X(1,inliers);X(4,inliers)+o], [X(2,inliers);X(5,inliers)], 'color', 'g');
    
    if (plotOutliers)

        plot(a.X1(1, outliers),a.X1(2, outliers),'or', 'MarkerSize',4, 'linewidth', 1.2); 
        plot(a.X2(1, outliers)+o,a.X2(2, outliers),'or', 'MarkerSize',4, 'linewidth', 1.2);    
        
        line([X(1,outliers);X(4,outliers)+o], [X(2,outliers);X(5,outliers)], 'color', 'r');
    end
    else
        plot(a.X1(1, inliers),a.X1(2, inliers),'or', 'MarkerSize',4, 'linewidth', 1.2); 
        plot(a.X2(1, inliers)+o,a.X2(2, inliers),'or', 'MarkerSize',4, 'linewidth', 1.2);    
        line([X(1,inliers);X(4,inliers)+o], [X(2,inliers);X(5,inliers)], 'color', 'r');
        plot(a.X1(1, outliers),a.X1(2, outliers),'or', 'MarkerSize',4, 'linewidth', 1.2); 
        plot(a.X2(1, outliers)+o,a.X2(2, outliers),'or', 'MarkerSize',4, 'linewidth', 1.2);    
        
        line([X(1,outliers);X(4,outliers)+o], [X(2,outliers);X(5,outliers)], 'color', 'r');
    end

    axis image off;
    drawnow;
    
end
