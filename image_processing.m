clear all;

% Processing camera setup file
digit_end_position_X = 0;
digit_end_position_Y = 0;
digit_end_position_Z = 0;
digit_start_position_X = 0;
digit_start_position_Y = 0;
digit_start_position_Z = 0;
digit_solution_X = 0;
digit_solution_Y = 0;
digit_solution_Z = 0;
digit_interval_Z = 0;
fid = fopen('TIFF/wt.txt');
tline = fgetl(fid);
while ischar(tline)
    matches = strfind(tline, 'AxisCode="X"');
    num = length(matches);
    if num > 0
        %fprintf(1,'%s\n',tline);
        tline = fgetl(fid);
        tline = fgetl(fid);
        tline = fgetl(fid);
        tline = fgetl(fid);
        % End position
        tline = fgetl(fid);
        fprintf(1,'%s\n',tline);
        digit_end_position_X = sscanf(tline, '%*[EndPosition=]%f');
        digit_end_position_X
        tline = fgetl(fid);
        digit_solution_X = sscanf(tline, '%*[GUI MaxSize=]%f');
        digit_solution_X
        tline = fgetl(fid);
        tline = fgetl(fid);
        tline = fgetl(fid);
        fprintf(1,'%s\n',tline);
        digit_start_position_X = sscanf(tline, '%*[StartPosition=]%f');
        digit_start_position_X
    end
    tline = fgetl(fid);
end
fclose(fid);

fid = fopen('TIFF/wt.txt');
tline = fgetl(fid);
while ischar(tline)
    matches = strfind(tline, 'AxisCode="Y"');
    num = length(matches);
    if num > 0
        %fprintf(1,'%s\n',tline);
        tline = fgetl(fid);
        tline = fgetl(fid);
        tline = fgetl(fid);
        tline = fgetl(fid);
        % End position
        tline = fgetl(fid);
        fprintf(1,'%s\n',tline);
        digit_end_position_Y = sscanf(tline, '%*[EndPosition=]%f');
        digit_end_position_Y
        tline = fgetl(fid);
        digit_solution_Y = sscanf(tline, '%*[GUI MaxSize=]%f');
        digit_solution_Y
        tline = fgetl(fid);
        tline = fgetl(fid);
        tline = fgetl(fid);
        fprintf(1,'%s\n',tline);
        digit_start_position_Y = sscanf(tline, '%*[StartPosition=]%f');
        digit_start_position_Y
    end
    tline = fgetl(fid);
end
fclose(fid);

fid = fopen('TIFF/wt.txt');
tline = fgetl(fid);
while ischar(tline)
    matches = strfind(tline, 'AxisCode="Z"');
    num = length(matches);
    if num > 0
        %fprintf(1,'%s\n',tline);
        tline = fgetl(fid);
        tline = fgetl(fid);
        tline = fgetl(fid);
        % End position
        tline = fgetl(fid);
        fprintf(1,'%s\n',tline);
        digit_end_position_Z = sscanf(tline, '%*[EndPosition=]%f');
        digit_end_position_Z
        tline = fgetl(fid);
        fprintf(1,'%s\n',tline);
        digit_solution_Z = sscanf(tline, '%*[GUI MaxSize=]%f');
        digit_solution_Z
        tline = fgetl(fid);
        fprintf(1,'%s\n',tline);
        digit_interval_Z = sscanf(tline, '%*[Interval=]%f');
        digit_interval_Z
        tline = fgetl(fid);
        tline = fgetl(fid);
        tline = fgetl(fid);
        tline = fgetl(fid);
        fprintf(1,'%s\n',tline);
        digit_start_position_Z = sscanf(tline, '%*[StartPosition=]%f');
        digit_start_position_Z
    end
    tline = fgetl(fid);
end
fclose(fid);


digitx = (digit_end_position_X - digit_start_position_X)/(digit_solution_X-1);
digity = (digit_end_position_Y - digit_start_position_Y)/(digit_solution_Y-1);
digitz = (digit_end_position_Z - digit_start_position_Z)/(digit_solution_Z-1);
digitx0 = digit_start_position_X;
digity0 = digit_start_position_Y;
digitz0 = digit_start_position_Z;
digitz = digitz/1000;
digitz0 = digitz0/1000;

final_L3D_in = [];

for T = 44:44
    final_L3D_in = [];
    
    filename = sprintf('%s%03d%s', 'out/point_outer_', T, '.csv');
    fileID = fopen(filename,'w');
    fprintf(fileID,'%s, %s, %s, %s\n', 'x coord', 'y coord', 'z coord', 'scalar');
    fclose(fileID);
    
    filename = sprintf('%s%03d%s', 'out/point_inner_', T, '.csv');
    fileID = fopen(filename,'w');
    fprintf(fileID,'%s, %s, %s, %s\n', 'x coord', 'y coord', 'z coord', 'scalar');
    fclose(fileID);
    
    for Z = 25:25
        filename = sprintf('%s%02d%s%03d%s', 'TIFF/wt_C001Z0', Z, 'T', T, '.tif');
        BW1 = imread(filename);
        
        BW2 = imadjust(BW1);
        A = double(BW2);
        B = zeros(512);
        stick_L = 35;
        stick_thres1 = 2800;
        stick_thres2 = 3000;
        %stick_thres = .03*max(max(A));
        %A(1:50, 1:50) = stick_thres;
        
        sigma=1.5;     % scale parameter in Gaussian kernel
        G=fspecial('gaussian',15,sigma);
        A=conv2(A,G,'same');  % smooth image by Gaussiin convolution
        
        %imshow(A');
        %plot(A)
        hold on
        for n = 1:1000000
            rline = randi(512,2);
            k = (rline(2,2) - rline(2,1))/(rline(1,2) - rline(1,1));
            rline(2,1) = rline(1,1) + stick_L/(sqrt(1+k*k));
            rline(2,2) = rline(1,2) + k*stick_L/(sqrt(1+k*k));
            %rline;
            use_flag1 = 0;
            x0 = round(rline(1,1));
            y0 = round(rline(1,2));
            if ( x0 > 0 && x0 < 512 )
                if ( y0 > 0 && y0 < 512 )
                    if ( A(x0, y0) > stick_thres1 )
                        use_flag1 = 1;
                    end
                end
            end
            use_flag2 = 0;
            x1 = round(rline(2,1));
            y1 = round(rline(2,2));
            if ( x1 > 0 && x1 < 512 )
                if ( y1 > 0 && y1 < 512 )
                    if ( A(x1, y1) > stick_thres1 )
                        use_flag2 = 1;
                    end
                end
            end
            if ( use_flag1 * use_flag2 == 1 )
                if ( (A(x1, y1) +  A(x0, y0))/2 > stick_thres2 )
                    %plot([x0,x1], [y0,y1])
                    for xpix = x0:x1
                        ypix = round(y0 + (xpix - x0)*k);
                        if (ypix >0 && ypix < 512)
                            B(xpix, ypix) = 1;
                        end
                    end
                end
            end
        end
        % select four point that most likely to be in the main body
        c = [0,0,0,0];
        r = [0,0,0,0];
        xc = 128;
        yc = 128;
        for i = 1:128
            if ( B(xc + i, yc+i) == 1)
                c(1) = xc+i;
                r(1) = yc+i;
                break;
            end
        end
        xc = 128+256;
        yc = 128;
        for i = 1:128
            if ( B(xc - i, yc+i) == 1)
                c(2) = xc-i;
                r(2) = yc+i;
                break;
            end
        end
        xc = 128;
        yc = 128+256;
        for i = 1:128
            if ( B(xc + i, yc-i) == 1)
                c(3) = xc+i;
                r(3) = yc-i;
                break;
            end
        end
        xc = 128+256;
        yc = 128+256;
        for i = 1:128
            if ( B(xc - i, yc-i) == 1)
                c(4) = xc-i;
                r(4) = yc-i;
                break;
            end
        end
        B2 = bwselect(B,r, c, 8);
        
        clf;
        imshow(B);
        outname = sprintf('%s%02d%s%03d%s', 'out/L', Z, 'T', T, '.tif');
        print('-dtiff',outname)
        
        clf;
        imshow(B2);
        outname = sprintf('%s%02d%s%03d%s', 'out/S', Z, 'T', T, '.tif');
        print('-dtiff',outname)
        
        B2 = double(B2);
        B3 = conv2(B2,G,'same');
        dx = 1;
        dy = 1;
        b = 0.3*ones(size(B2));
        B3 = evolve2D(B3,dx,dy,0.5,1000,[],[],0,[],0,[],[],1,b);
        
%         % do something to the boundary so that we dont get broken contours
%         B3(B3(1:2, :) > 0.5) = 0.5;
%         B3(B3(511:512, :) > 0.5) = 0.5;
%         B3(B3(:, 1:2) > 0.5)= 0.5;
%         B3(B3(:, 511:512) > 0.5) = 0.5;
        L = contour(B3, [0.5,0.5]);
        B3T = B3';
        
        %L3D = L;
        
        hold on
        filename = sprintf('%s%03d%s', 'out/point_outer_', T, '.csv');
        fileID1 = fopen(filename,'a');
        filename = sprintf('%s%03d%s', 'out/point_inner_', T, '.csv');
        fileID2 = fopen(filename,'a');
        
        clf;
        hold on;
        loop_index=1;
        while ( loop_index<length(L) )
            %L3D(1, loop_index) = digitz0+(Z-0.5)*digitz;
            loop_length = L(2,loop_index);
            loop_index = loop_index + 1;
            LL = L(:, loop_index:loop_index + loop_length - 1);
            
            % finding a point on the contour to check in-or-out
            find_check_point = 0;
            % Searching left boundary
            p1 = find(LL(1,:) == min(LL(1,:)), 1);
            x1 = LL(1,p1);
            y1 = LL(2,p1);
            [LL_sort, IX] = sort(LL(1,:));
            p2 = 2;
            x2 = LL(1,IX(p2));
            y2 = LL(2,IX(p2));
            while (y2 > y1)
                p2 = p2 + 1;
                x2 = LL(1,IX(p2));
                y2 = LL(2,IX(p2));
            end
            dx = -(y2-y1);
            dy = x2-x1;
            x2 = round(x1 + 10*dx);
            y2 = round(y1 + 10*dy);
            if ( x2 > 0 && x2 < 512 )
                if (y2 > 0 && y2 < 512)
                    if ( dx^2 + dy^2 < 4 )
                    find_check_point = 1;
                    plot(x1,y1,'*g','markersize', 6)
                    plot(x2,y2,'*r','markersize', 6)
                    end
                end
            end
            
            if (find_check_point ==0)
                %searching down
                p1 = find(LL(2,:) == min(LL(2,:)), 1);
                x1 = LL(1,p1);
                y1 = LL(2,p1);
                [LL_sort, IX] = sort(LL(2,:));
                p2 = 2;
                x2 = LL(1,IX(p2));
                y2 = LL(2,IX(p2));
                while (x2 < x1)
                    p2 = p2 + 1;
                    x2 = LL(1,IX(p2));
                    y2 = LL(2,IX(p2));
                end
                dx = -(y2-y1);
                dy = x2-x1;
                x2 = round(x1 + 10*dx);
                y2 = round(y1 + 10*dy);
                if ( x2 > 0 && x2 < 512 && y2 > 0 && y2 < 512)
                    if ( dx^2 + dy^2 < 4 )
                    find_check_point = 1;
                    plot(x1,y1,'*g','markersize', 6)
                    plot(x2,y2,'*r','markersize', 6)
                    end
                end
            end
            
            if (find_check_point ==0)
                % Searching right
                p1 = find(LL(1,:) == max(LL(1,:)), 1);
                x1 = LL(1,p1);
                y1 = LL(2,p1);
                [LL_sort, IX] = sort(LL(1,:));
                p2 = length(IX) - 1;
                x2 = LL(1,IX(p2));
                y2 = LL(2,IX(p2));
                while (y2 < y1)
                    p2 = p2 - 1;
                    x2 = LL(1,IX(p2));
                    y2 = LL(2,IX(p2));
                end
                dx = -(y2-y1);
                dy = x2-x1;
                x2 = round(x1 + 10*dx);
                y2 = round(y1 + 10*dy);
                if ( x2 > 0 && x2 < 512 && y2 > 0 && y2 < 512)
                    if ( dx^2 + dy^2 < 4 )
                    find_check_point = 1;
                    plot(x1,y1,'*g','markersize', 6)
                    plot(x2,y2,'*r','markersize', 6)
                    end
                end
            end
            
            if (find_check_point ==0)
                % Searching top
                p1 = find(LL(2,:) == max(LL(2,:)), 1);
                x1 = LL(1,p1);
                y1 = LL(2,p1);
                [LL_sort, IX] = sort(LL(2,:));
                p2 = length(IX) - 1;
                x2 = LL(1,IX(p2));
                y2 = LL(2,IX(p2));
                while (x2 > x1)
                    p2 = p2 - 1;
                    x2 = LL(1,IX(p2));
                    y2 = LL(2,IX(p2));
                end
                dx = -(y2-y1);
                dy = x2-x1;
                x2 = round(x1 + 10*dx);
                y2 = round(y1 + 10*dy);
                if ( x2 > 0 && x2 < 512 && y2 > 0 && y2 < 512)
                    if ( dx^2 + dy^2 < 4 )
                    find_check_point = 1;
                    plot(x1,y1,'*g','markersize', 6)
                    plot(x2,y2,'*r','markersize', 6)
                    end
                end
            end
            % transform
            if (find_check_point == 1)
                if ( B3T(x2, y2) > 0.5 )
                    loop_color = 1;
                else
                    loop_color = 0;
                end
            else
                loop_color = 1;
            end
            if (loop_color == 1)
                fileID = fileID1;
            else
                fileID = fileID2;
            end
            for i=1:loop_length
                fprintf(fileID,'%6.2f, %6.2f, %6.2f, %6.2f\n', digitx0+(L(1,loop_index)-0.5)*digitx, digity0+(L(2,loop_index)-0.5)*digity, digitz0+(Z-0.5)*digitz, loop_color);
                loop_index = loop_index + 1;
            end
            if ( loop_color == 0 )
                if (length(final_L3D_in) == 0)
                    final_L3D_in = [[digitz0+(Z-0.5)*digitz, loop_length]' LL];
                else
                    final_L3D_in = [final_L3D_in [[digitz0+(Z-0.5)*digitz, loop_length]' LL]];
                end
            end
            
            %             % finding a point on the contour to check in-or-out
            %             find_check_point = 0;
            %             % Searching left boundary
            %             p1 = find(LL(1,:) == min(LL(1,:)), 1);
            %             x1 = LL(1,p1);
            %             y1 = LL(2,p1);
            %             [LL_sort, IX] = sort(LL(1,:));
            %             p2 = 2;
            %             x2 = LL(1,IX(p2));
            %             y2 = LL(2,IX(p2));
            %             while (y2 > y1)
            %                 p2 = p2 + 1;
            %                 x2 = LL(1,IX(p2));
            %                 y2 = LL(2,IX(p2));
            %             end
            %             dx = -(y2-y1);
            %             dy = x2-x1;
            %             x2 = round(x1 + 10*dx);
            %             y2 = round(y1 + 10*dy);
            %             if ( x2 > 0 && x2 < 512 )
            %                 if (y2 > 0 && y2 < 512)
            %                     find_check_point = 1;
            %                     plot(x1,y1,'*g','markersize', 6)
            %                     plot(x2,y2,'*r','markersize', 6)
            %            % transform
            %             if ( B3T(x2, y2) > 0.5 )
            %                 loop_color = 1;
            %             else
            %                 loop_color = 0;
            %             end
            %                 end
            %             end
            %             if (find_check_point == 0)
            %                 loop_color = 1;
            %             end
            %             if (loop_color == 1)
            %                 fileID = fileID1;
            %             else
            %                 fileID = fileID2;
            %             end
            %             for i=1:loop_length
            %                 fprintf(fileID,'%6.2f, %6.2f, %6.2f, %6.2f\n', digitx0+(L(1,loop_index)-0.5)*digitx, digity0+(L(2,loop_index)-0.5)*digity, digitz0+(Z-0.5)*digitz, loop_color);
            %                 loop_index = loop_index + 1;
            %             end
            %             if ( loop_color == 0 )
            %                 if (length(final_L3D_in) == 0)
            %                     final_L3D_in = [[Z, loop_length]' LL];
            %                 else
            %                     final_L3D_in = [final_L3D_in [[Z, loop_length]' LL]];
            %                 end
            %             end
        end
        fclose(fileID1);
        fclose(fileID2);
        %pause(10);
        
        clf;
        imshow(BW2);
        hold on;
        for i = 1:length(L)
            plot(L(1,i),L(2,i), 'color', 'green');
        end
        
        outname = sprintf('%s%02d%s%03d%s', 'out/Z', Z, 'T', T, '.tif');
        print('-dtiff',outname)
        
        clf;
        %figure;
        imshow(B3);
        outname = sprintf('%s%02d%s%03d%s', 'out/B', Z, 'T', T, '.tif');
        print('-dtiff',outname)
        
        %
        %         % select four point that most likely to be in the main body
        %         c = [0,0,0,0];
        %         r = [0,0,0,0];
        %         xc = 128;
        %         yc = 128;
        %         for i = 1:128
        %             if ( B(xc + i, yc+i) == 1)
        %                 c(1) = xc+i;
        %                 r(1) = yc+i;
        %                 break;
        %             end
        %         end
        %         xc = 128+256;
        %         yc = 128;
        %         for i = 1:128
        %             if ( B(xc - i, yc+i) == 1)
        %                 c(2) = xc-i;
        %                 r(2) = yc+i;
        %                 break;
        %             end
        %         end
        %         xc = 128;
        %         yc = 128+256;
        %         for i = 1:128
        %             if ( B(xc + i, yc-i) == 1)
        %                 c(3) = xc+i;
        %                 r(3) = yc-i;
        %                 break;
        %             end
        %         end
        %         xc = 128+256;
        %         yc = 128+256;
        %         for i = 1:128
        %             if ( B(xc - i, yc-i) == 1)
        %                 c(4) = xc-i;
        %                 r(4) = yc-i;
        %                 break;
        %             end
        %         end
        %         B2 = bwselect(B,r, c, 8);
        %
        %         B2 = double(B2);
        %         B3 = conv2(B2,G,'same');
        %         dx = 1;
        %         dy = 1;
        %         b = 0.3*ones(size(B2));
        %         B3 = evolve2D(B3,dx,dy,0.5,1000,[],[],0,[],0,[],[],1,b);
        %         L = contour(B3, [0.5,0.5]);
        %         B3T = B3';
        %
        %         L3D = L;
        %
        %         hold on
        %         filename = sprintf('%s%03d%s', 'out/point_outer_', T, '.csv');
        %         fileID1 = fopen(filename,'a');
        %         filename = sprintf('%s%03d%s', 'out/point_inner_', T, '.csv');
        %         fileID2 = fopen(filename,'a');
        %
        %         loop_index=1;
        %         while ( loop_index<length(L) )
        %             L3D(1, loop_index) = Z;
        %             loop_length = L(2,loop_index);
        %             loop_index = loop_index + 1;
        %             LL = L(:, loop_index:loop_index + loop_length - 1);
        %             p1 = find(LL(1,:) == min(LL(1,:)));
        %             x1 = LL(1,p1);
        %             y1 = LL(2,p1);
        %             [LL_sort, IX] = sort(LL(1,:));
        %             p2 = 2;
        %             x2 = LL(1,IX(p2));
        %             y2 = LL(2,IX(p2));
        %             while (y2 > y1)
        %                 p2 = p2 + 1;
        %                 x2 = LL(1,IX(p2));
        %                 y2 = LL(2,IX(p2));
        %             end
        %             plot(x1,y1,'*g','markersize', 6)
        %             plot(x2,y2,'*r','markersize', 6)
        %             dx = -(y2-y1);
        %             dy = x2-x1;
        %             x2 = round(x1 + 10*dx);
        %             y2 = round(y1 + 10*dy);
        %             if ( B3T(x2, y2) > 0.5 )
        %                 loop_color = 1;
        %             else
        %                 loop_color = 0;
        %             end
        %             if (loop_color == 1)
        %                 fileID = fileID1;
        %             else
        %                 fileID = fileID2;
        %             end
        %             for i=1:loop_length
        %                 fprintf(fileID,'%6.2f, %6.2f, %6.2f, %6.2f\n', digitx0+(L(1,loop_index)-0.5)*digitx, digity0+(L(2,loop_index)-0.5)*digity, digitz0+(Z-0.5)*digitz, loop_color);
        %                 loop_index = loop_index + 1;
        %             end
        %             if ( loop_color == 0 )
        %                 if (length(final_L3D_in) == 0)
        %                     final_L3D_in = [[Z, loop_length]' LL];
        %                 else
        %                     final_L3D_in = [final_L3D_in [[Z, loop_length]' LL]];
        %                 end
        %             end
        %         end
        %         fclose(fileID1);
        %         fclose(fileID2);
        %
        %         imshow(BW2);
        %         hold on;
        %         for i = 1:length(L)
        %             plot(L(1,i),L(2,i), 'color', 'green');
        %         end
        %
        %         outname = sprintf('%s%02d%s%03d%s', 'out/Z', Z, 'T', T, '.tif');
        %         print('-dtiff',outname)
        %
        %         clf;
        %         %figure;
        %         imshow(B3);
        %         outname = sprintf('%s%02d%s%03d%s', 'out/B', Z, 'T', T, '.tif');
        %         print('-dtiff',outname)
        
        % % This method does not work!!!
        % % scan the outer box to find a boundary of the tube
        % scan_offset = 50;
        % scan_outer = 0.01;
        % scan_flag = 1;
        % scan_direction = -1;
        % for i = 1:512
        %     if ( B3T(scan_offset, i) > 0.01 )
        %         scan_flag = 0;
        %         break;
        %     end
        %     scan_direction = 1;
        % end
        % for i = 1:512
        %     if ( B3T(i, scan_offset) > 0.01 )
        %         scan_flag = 0;
        %         break;
        %     end
        %     scan_direction = 2;
        % end
        % for i = 1:512
        %     if ( B3T(512-scan_offset, i) > 0.01 )
        %         scan_flag = 0;
        %         break;
        %     end
        %     scan_direction = 3;
        % end
        % for i = 1:512
        %     if ( B3T(i, 512-scan_offset) > 0.01 )
        %         scan_flag = 0;
        %         break;
        %     end
        %     scan_direction = 4;
        % end
        %
        % imshow(B3);
        % hold on
        % switch(scan_direction)
        %     case 1
        %
        %     case 2
        %
        %     case 3
        %
        %     case 4
        %         for i = 1:512
        %             scan_x = i;
        %             level_mark = -1;
        %             for j = 1:512
        %                 scan_y = 513 - j;
        %                 if ( level_mark*(B3T(scan_x, scan_y) - 0.5) < 0 )
        %                     if (level_mark > 0)
        %                         plot(scan_x, scan_y, 'r');
        %                     else
        %                         plot(scan_x, scan_y, 'g');
        %                     end
        %                     level_mark = -1*level_mark;
        %                 end
        %             end
        %         end
        %     otherwise
        %         disp('wrong scan_direction...stop');
        %         return
        % end
        
    end
    filename = sprintf('%s%03d%s', 'out/LayeredContour', T, 'txt');
    save(filename, 'final_L3D_in', '-ASCII')
end
