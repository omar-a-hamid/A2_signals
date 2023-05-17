% Read the input image
image = imread('Ocean.bmp');

% Define the parameters
A1 = 30;
B1 = 20;
C1 = 180;
D1 = 230;

A2 = 70;
B2 = 20;
C2 = 140;
D2 = 240;
image_path = 'Ocean.bmp';



A = A2;
B = B2;
C = C2;
D = D2;

ip_axis = 0:255;
op_1 = ip_axis(ip_axis<=A).*B./A;
% op_2 = (ip_axis((A<ip_axis<C)&&(ip_axis<C))).*((D-B)./(C-A))+B;(img(img >= r1 & img <= r2) 
op_2 = ((ip_axis(A<=ip_axis & ip_axis<=C)-A).*((D-B)./(C-A)))+B;
op_3 = (ip_axis(C<=ip_axis)-C).*((255-D)./(255-C))+D;
figure;
hold on 
grid on
plot(ip_axis(1:length(op_1)),op_1)
plot(ip_axis(A<=ip_axis & ip_axis<=C),op_2)
plot(ip_axis(C<=ip_axis),op_3)
%%
% Transformation parameters for Ocean_1.bmp

transformed_image1=contrast_fn(image_path, A1, B1, C1, D1);

% Transformation parameters for Ocean_2.bmp

transformed_image2=contrast_fn(image_path, A2, B2, C2, D2);
transformed_image3=contrast_fn(image_path, 20, 20, 200, 200); %linear should return original 
transformed_image4=contrast_fn(image_path, 125, 0, 125, 255); %binary image
imwrite(transformed_image1, 'Ocean_1.bmp');
imwrite(transformed_image2, 'Ocean_2.bmp');
imwrite(transformed_image3, 'Ocean_3.bmp');
imwrite(transformed_image4, 'Ocean_4.bmp');



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         functions       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
function transformed_image=contrast_fn(image_path, r1, s1, r2, s2)
    % Load the grayscale image
    img = imread(image_path);

    % Normalize image intensity values
    img = im2double(img);

    % Contrast stretching transformation
    out = zeros(size(img));

    out(img < r1) = s1 / r1 * img(img < r1);
    out(img >= r1 & img <= r2) = ((s2 - s1) / (r2 - r1)) * (img(img >= r1 & img <= r2) - r1) + s1;
    out(img > r2) = ((1 - s2) / (1 - r2)) * (img(img > r2) - r2) + s2;

    % Convert back to original intensity range
    transformed_image = im2uint8(out);

    % Save the transformed image
    [~, name, ~] = fileparts(image_path);
    transformed_path = [name '_transformed.bmp'];
    imwrite(transformed_image, transformed_path);

    fprintf('Transformed image saved as: %s\n', transformed_path);
end
