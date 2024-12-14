% Load a colorful image
imagePath = 'example_image.jpg'; % Replace with your image file path
img = imread(imagePath);

% Normalize the image RGB values to [0, 1]
imgNormalized = double(img) / 255;

% Define the custom colormap for heights (20nm to 220nm)
% Example: Heights (in nm) and their corresponding RGB values
numHeights = 200;
minHeight = 20; % Minimum height in nm
maxHeight = 220; % Maximum height in nm
heights = linspace(minHeight, maxHeight, numHeights); % Heights from 20nm to 220nm

% Mock example of RGB values for heights, replace with real interference color data
customColormap = jet(numHeights); % Replace with your experimental RGB data

% Ensure the colormap is in the range [0, 1]
customColormap = double(customColormap);

% Convert RGB image to indexed image based on the custom colormap
[indexedImage, ~] = rgb2ind(imgNormalized, customColormap);

% Map the indices to heights
heightMap = heights(indexedImage + 1); % Add 1 because MATLAB indices start at 1

% Smooth the height map
smoothedHeightMap = imgaussfilt(heightMap, 2); % Gaussian smoothing with sigma = 2

% Visualize the height map
figure;
surf(flipud(smoothedHeightMap), 'EdgeColor', 'none'); % Create a surface plot
colormap(jet); % Example colormap for visualization; adjust if needed
colorbar; % Add a colorbar to visualize height range
axis tight; % Adjust axis
view(2); % Top-down 2D view

% Add titles and labels
title('Height Map from Interference Colors (Using rgb2ind)');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Height (nm)');

% Enhance visualization
shading interp; % Smooth shading