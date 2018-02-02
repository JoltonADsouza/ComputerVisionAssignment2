clc
close all;
clear all;

GroundTruth = imread('C:\Users\Jolton\Desktop\Files\Homeworks\Semester3\3D Computer Vision\Homework2\disp2.pgm');
figure
% GroundTruth = int8(GroundTruth);
imshow(GroundTruth);

I_L = imread('C:\Users\Jolton\Desktop\Files\Homeworks\Semester3\3D Computer Vision\Homework2\teddyL.pgm');
figure
% I_L = int8(I_L);
imshow(I_L);

I_R = imread('C:\Users\Jolton\Desktop\Files\Homeworks\Semester3\3D Computer Vision\Homework2\teddyR.pgm');
figure
% I_R = int8(I_R);
imshow(I_R);

I_L_ranked = zeros(size(I_L,1),size(I_L,2));
I_R_ranked = zeros(size(I_R,1),size(I_R,2));

I_L = padarray(I_L,[2 2],0,'both');

I_R = padarray(I_R,[2 2],0,'both');

rnk1 = 0;
for i = 3:size(I_L,1)-2
    for j = 3:size(I_L,2)-2
        for k = i-2:i+2
            for l = j-2:j+2
                if I_L(k,l) < I_L(i,j)
                    rnk1 = rnk1+1;
                end
            end
        end
        I_L_ranked(i-2,j-2) = rnk1;
        rnk1 = 0;
    end
end

I_L_ranked = padarray(I_L_ranked,[7 7],0,'both');


rnk2 = 0;
for i = 3:size(I_R,1)-2
    for j = 3:size(I_R,2)-2
        for k = i-2:i+2
            for l = j-2:j+2
                if I_R(k,l) < I_R(i,j)
                    rnk2 = rnk2+1;
                end
            end
        end
        I_R_ranked(i-2,j-2) = rnk2;
        rnk2 = 0;
    end
end
I_R_ranked = padarray(I_R_ranked,[7 7],0,'both');



minimum = 255; brk = 0; c = 1; d = 1; sum = 0; o =1; p = 1; q = 1; SAD_array = zeros(1,63); disparity_map = zeros(1,2); SAD_min = zeros(1,2);
for i = 1:size(I_R_ranked,1)-14
    c = i-1;
       for j = 1:size(I_R_ranked,2)-14
        d = j-1;
        x = j+62;
        if x > size(I_R_ranked,2)-14
           x = size(I_R_ranked,2)-14;
           SAD_array = zeros(1,63-((j+62)-(size(I_R_ranked,2)-7)));
        end
        for k = i:size(I_L_ranked,1)-14
            for l = j:x
                for m = k:k+14
                    c = c+1;
                    for n = l:l+14
                        d = d+1;
                         difference = int16(I_R_ranked(c,d)) - int16(I_L_ranked(m,n));                           
                         sum = sum + abs(difference);
                    end
                        d = j-1;
                end

                    SAD_array(o) = sum;
                    o = o+1;
                    sum = 0;
%                     c1 = c
                    c = i-1;
            end
                
                SAD_min(q,p) = min(SAD_array); 
                [y,SAD_min_index] = find(SAD_array == min(SAD_array));
                disparity_map(q,p) = (j+(min(SAD_min_index)-1)) - p;

                p = p+1;
                SAD_array = zeros(1,63);
                o = 1;
                if j == size(I_R_ranked,2)-14
                    brk = brk+1;
                    break;
                else
                    break;
                end
        end
            if brk == 1
                p = 1;
                q = q+1;
                brk = 0;
                break;
            end
      end
end

figure;
imshow(uint8(disparity_map*4));

%% Error Rate
bad_pixel_count = 0;
GroundTruth = round(GroundTruth/4);
for i = 1:size(GroundTruth,1)
    for j = 1:size(GroundTruth,2)
        if abs(disparity_map(i,j) - GroundTruth(i,j)) > 1
            bad_pixel_count = bad_pixel_count + 1;
        end
    end
end
error_percent = (bad_pixel_count/   (size(GroundTruth,1)*size(GroundTruth,2)))*100