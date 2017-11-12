hold on;
%for T = 2:47
T = 2;
filename = sprintf('%s%03d%s', 'out/point_inner_', T, '.csv');
    Cin = csvread(filename,1);
    TPnum = length(Cin(:,1));
    
    x0 = 521/2;
    y0 = 512/2;
    R = sqrt(x0^2 + y0^2);
    
    dspan_min = 512*15;
    for ktheta = 0.1:0.1:2*pi
        xb = x0 + R*cos(ktheta);
        yb = y0 + R*sin(ktheta);
        kline = -1/tan(ktheta);
        
        dmax = 0;
        dmin = 512*1.5;
        for i = 1: TPnum
            xp = Cin(i, 1);
            yp = Cin(i, 2);
            d = abs(yp - yb - kline*(xp - xb))/sqrt(kline^2 + 1);
            if ( d > dmax)
                dmax = d;
                imax = i;
            end
            if ( d < dmin)
                dmin = d;
                imin = i;
            end
        end
        
        dspan = dmax - dmin;
        if ( dspan < dspan_min )
            dspan_min = dspan;
            ispan1 = imax;
            ispan2 = imin;
            kspan = kline;
            theta_cut = ktheta + pi/2;
        end
        %plot(Cin(imax, 1), Cin(imax, 2), 'or');
        %plot(Cin(imin, 1), Cin(imin, 2), 'or');
    end
    kcut = -1/kspan;
    % save cut slope
    filename = sprintf('%s%03d%s', 'out/cutslope_', T, '.csv');
    fileID = fopen(filename,'w');
    fprintf(fileID,'%f\n', kcut);
    fclose(fileID);
    
    
    filename = sprintf('%s%03d%s', 'out/LayeredContour', T, 'txt');
    %load(filename, 'final_L3D_in');
    final_L3D_in = load(filename);
    Lrun = linspace(-R,R,300);
    Xrun = linspace(0,512);
    Lnew = final_L3D_in;
    %Lnew = filename;
    
    %maxial = [];
    
    filename = sprintf('%s%03d%s', 'out/central_', T, '.csv');
    fileID = fopen(filename,'w');
    fprintf(fileID,'%s, %s, %s, %s\n', 'x coord', 'y coord', 'z coord', 'scalar');
    
    for i = 1:length(Lrun)
        %for i = 1:1
        x1 = x0 + Lrun(i)*cos(theta_cut);
        y1 = y0 + Lrun(i)*sin(theta_cut);
        Yrun = y1 + kcut*(Xrun-x1);
        
        loop_index=1;
        dist_cut_square_max = 0;
        mx = 0;
        my = 0;
        mz = 0;
        nlayer = 0;
        while ( loop_index < length(Lnew) )
            loop_length = Lnew(2,loop_index);
            loop_index = loop_index + 1;
            LL = Lnew(:, loop_index:loop_index + loop_length - 1);
            
            plot(LL(1,:), LL(2,:))
            %plot(Cupper(:,1), Cupper(:,2),'r')
            
            PSect = InterX([Xrun;Yrun],[LL(1,:);LL(2,:)]);
            if ( length(PSect) > 2 )
                disp('warning: more than 3 intersection points are found!');
            end
            
            if (length(PSect(1,:)) == 2)
                %dist_cut_square = (PSect(1,2) - PSect(1,1))^2 + (PSect(2,2) - PSect(2,1))^2;
                %if ( dist_cut_square > dist_cut_square_max )
                %dist_cut_square_max = dist_cut_square;
                mx = mx + (PSect(1,2) + PSect(1,1))/2;
                my = my + (PSect(2,2) + PSect(2,1))/2;
                mz = mz + Lnew(1, loop_index-1);
                nlayer = nlayer + 1;
                %end
            end
            %plot(PSect(1,:), PSect(2,:),'');
            loop_index = loop_index + loop_length;
        end
        if (nlayer > 0)
            plot(mx/nlayer, my/nlayer, 'm*')
            fprintf(fileID,'%6.2f, %6.2f, %6.2f\n', mx/nlayer, my/nlayer, mz/nlayer);
        end
        plot(Xrun,Yrun,'-.k');
    end
    
    fclose(fileID);
    
    % for Z = 7:7
    %     zlower = Z;
    %     zupper = Z+1;
    %     Clower = Cin(Cin(:,3)==zlower, :);
    %     Cupper = Cin(Cin(:,3)==zupper, :);
    %
    %     plot(Clower(:,1), Clower(:,2))
    %     %plot(Cupper(:,1), Cupper(:,2),'r')
    %
    %     for i = 1:length(Lrun)
    %         x1 = x0 + Lrun(i)*cos(theta_cut);
    %         y1 = y0 + Lrun(i)*sin(theta_cut);
    %         Yrun = y1 + kcut*(Xrun-x1);
    %         PSect = InterX([Xrun;Yrun],[Clower(:,1)';Clower(:,2)']);
    %         if ( length(PSect) > 2 )
    %             disp('warning:');
    %         end
    %         plot(PSect(1,:), PSect(2,:),'m*');
    %         plot(Xrun,Yrun,'-.k');
    %     end
    %
    % %     [vx,vy]=voronoi([Clower(1:10:length(Clower(:,1)),1); Cupper(1:10:length(Cupper(:,1)),1)],[Clower(1:10:length(Clower(:,1)),2); Cupper(1:10:length(Cupper(:,1)),2)]);
    % %     plot(vx, vy,'-*r');
    % %     xlim([0, 512]);
    % %     ylim([0, 512]);
    %
    %
    %     %xy = [0 0; 60 0; 100 20; 0 80];
    % %     xy1 = Clower(1:50:length(Clower(:,1)),1:2);
    % %     xy2 = Cupper(1:50:length(Cupper(:,1)),1:2);
    % %     plot(xy1(:,1), xy1(:,2));
    % %     plot(xy2(:,1), xy2(:,2), 'r');
    % %     pdeGeom = [geomDataFromPolygon(xy1) geomDataFromPolygon(xy2)];
    % %     [medialCurves,pdeGeomPadded,box] = computeMedialAxis(pdeGeom);
    % %     plotMedialCurves(medialCurves);
    %     %for i = 1:length(Clower)
    % end
    
    
    plot(Cin(ispan1, 1), Cin(ispan1, 2), 'sg');
    plot(Cin(ispan2, 1), Cin(ispan2, 2), 'sg');
    
    x = 1:512;
    line1 = Cin(ispan1,2) + kspan*(x - Cin(ispan1,1));
    line2 = Cin(ispan2,2) + kspan*(x - Cin(ispan2,1));
    plot(x, line1, 'k')
    plot(x, line2, 'k')
    xlim([0,512]);
    ylim([0,512]);
    T
%end