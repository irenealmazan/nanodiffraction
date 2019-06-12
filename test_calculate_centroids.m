function centroid_struct = test_calculate_centroids_fullRock_Curve(dataout,fly2Dmaps)
    % this funtion computes the centroid for each pixel of a rocking curve
    % (over many angles)

            tempim = zeros(size(dataout,1),size(dataout,2));
            tempxcen = zeros(size(dataout,1),size(dataout,2));
            tempycen = zeros(size(dataout,1),size(dataout,2));
            tempthcen = zeros(size(dataout,1),size(dataout,2));
            
            for kk = 1:size(dataout,1)
                for ll = 1:size(dataout,2)
                    tempim(kk,ll) = fly2Dmaps.ii(kk).jj(ll).SumInt; %sum_angle sum(sum(ccd))
                    imgin =fly2Dmaps.ii(kk).jj(ll).im; %sum_angle ccd_angle -> 512 x 512 image
                    imgin_theta = fly2Dmaps.ii(kk).jj(ll).intensity; % sum(sum(ccd))_angle -> array with number of angles in rock curve entries
                    
                    line1=sum(imgin,1);  % vertical
                    line2=sum(imgin,2);  % horizontal
                    for mm=1:size(line1,2)
                        tempycen(kk,ll)=tempycen(kk,ll)+mm*line1(mm)/tempim(kk,ll) ;
                    end
                    for mm=1:size(line2,1)
                        tempxcen(kk,ll)=tempxcen(kk,ll)+mm*line2(mm)/tempim(kk,ll);
                    end
                    
                    for tt = 1:size(imgin_theta)
                       tempthcen(kk,ll) = tempthcen(kk,ll) + tt*imgin_theta/tempim(kk,ll);
                    end
                end
            end
            
            centroid_struct.Xcentroids = tempxcen; % pixel units
            centroid_struct.Ycentroids = tempycen;
            centroid_struct.Thcentroids = tempthcen;


            centroid_struct.imapx = dataout(1,:,3);
            centroid_struct.imapy = dataout(:,1,2);
            
end