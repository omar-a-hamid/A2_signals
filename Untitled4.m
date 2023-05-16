%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         part 2       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Load and convert the image to grayscale

clear
close all
clc

reduction_factor =5;
scan_range=2;
image = imread('coins.jpg');
%img_bw = rgb2gray(image);
%change image to gray scale
img_bw = 0.2989 * image(:,:,1) + 0.5870 * image(:,:,2) + 0.1140 * image(:,:,3);

binary_image = zeros(size(img_bw));
% turn to black and white for value less than threshhol value

binary_image(img_bw>78) = 255;

figure;imshow(binary_image)
[height, width] = size(binary_image);
label_counter=0;
current_label =0;
label = containers.Map;

objectx =[];
objecty =[];

boundryx =[];
boundryy =[];

% label('1,7')=1;
% label('2,8')=2;
%
% isKey(label,('1,7'))
% label('1,7')
mod_binary_reduction_image = ones(size(binary_image));
mod_binary_reduction_image = mod_binary_reduction_image*255;

% for i = 1:reduction_factor
%     for j = 1:reduction_factor
%         mod_binary_reduction_image(i,:)=0;
%         mod_binary_reduction_image(height-i,:)=0;
%         mod_binary_reduction_image(:,width-j)=0;
%         mod_binary_reduction_image(:,j)=0;
%         
%     end
% end

constant_i =0;
constant_j =0;
figure;imshow(mod_binary_reduction_image)

%erousion to remove any unwanted noise
for i = 0+1:height
    for j = 0+1:width
        if binary_image(i,j)==0
            
            while(i-reduction_factor+constant_i)<1
                constant_i=constant_i+1;
            end
            while(j-reduction_factor+constant_j)<1
                constant_j=constant_j+1;
            end
            mod_binary_reduction_image(i-reduction_factor+constant_i:i+reduction_factor,j-reduction_factor+constant_j:j+reduction_factor)=0;
%              mod_binary_image(i,j)=255;
            constant_i =0;
            constant_j =0;

        end
    end
end
figure;imshow(mod_binary_reduction_image)
% smoothed = conv(mod_binary_reduction_image((1:end-1),:),mod_binary_reduction_image((2:end),:),'same');
% figure;imshow(smoothed)

mod_binary_image = zeros(size(binary_image));
%dialation to remove the effect of errousion on wanted objects
for i = 1:height
    for j = 1:width
        if mod_binary_reduction_image(i,j)==255
            
            mod_binary_image(i-reduction_factor:i+reduction_factor,j-reduction_factor:j+reduction_factor)=255;
% mod_binary_image(i-reduction_factor:i+reduction_factor,j)=255;
%              mod_binary_image(i,j)=255;
        end
    end
end

figure;imshow(mod_binary_image)
binary_image = mod_binary_image;
% binary_image = mod_binary_reduction_image;

%check for 8-adjancy
for i = 1:height
    for j = 1:width
        
        if binary_image(i,j) == 0
            continue
        end
        %if the pixel is white check its nighbours
        if binary_image(i,j)==255   
            current_label =0;
            for i_nieghbour = -scan_range:scan_range
                if current_label
                    break
                end
                for j_nieghbour = -scan_range:scan_range
                    %if any of its nighbours has label on it use the same label
                    if isKey(label,strcat(num2str(i+i_nieghbour),',',num2str(j+j_nieghbour)))
                        current_label = label(strcat(num2str(i+i_nieghbour),',',num2str(j+j_nieghbour)));
                        label(strcat(num2str(i),',',num2str(j))) = current_label;
                        %add coordinates of each maching label to an array
                        objectx(current_label, end+1) =  i;
                        objecty(current_label, end+1) = j;
                        break;
                        
                    end
                end
            end
            %if its nighbours dont have labels creat new label and increment the labebl counter
            if ~current_label
                label_counter = label_counter+1 ;
                label(strcat(num2str(i),',',num2str(j))) = label_counter;
                objectx(end+1,1)= [i];
                objecty(end+1,1) =  [j];
                
            end
        end
        
    end
end
% [object_number, elemntsX] = size(objectx);
% [elemntsY, object_number] = size(objecty);


% image_mod = zeros(size(image));
%%
centers =[];
size = [];
%take the average of the coordinates of every object to find its center 
for i = 1: label_counter

    current_x=find(objectx(i,:));
    current_y=find(objecty(i,:));
    centerx = round(mean(objectx(i,current_x(:))));
    size(end+1) = max(objectx(i,current_x(:))) - min(objectx(i,current_x(:)));
    
    centery = round(mean(objecty(i,current_y(:))));
    image(centerx-1:centerx+1,centery-1:centery+1,:)=250;
    
    centers(i,2) = centerx-5;%adjust for text box size
    centers(i,1) = centery-7;%adjust for text box size




end
%%normalize the size in refrence to the smallest object

size = size - min(size);
% size = size .* 1.0 ./ max(size)
small_coins = length(size(size<5)); %any value less than the threshhold is a small object 
value =0;
value = value + small_coins *50;
large_coins = label_counter-small_coins;
value = value + large_coins *100;
figure;imshow(image)


for ii=1:length(centers)
   text_str{ii} = ['Coin' num2str(ii)];
end
RGB = insertText(image,centers,text_str,'FontSize',5);
figure;imshow(RGB)
title("labeled Image")
imwrite(RGB,'coins labeled image.png')
fprintf(strcat('Number of coins: ',num2str(label_counter),'\n'))
fprintf(strcat('Value: ',num2str(value),'\n'))
