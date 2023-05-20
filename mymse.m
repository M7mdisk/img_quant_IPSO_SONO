function mses = mymse(colorPalettes, image)
    [palleteSize, numPalettes ]= size(colorPalettes);
    matrix = double(image);



    mat = reshape(matrix, [], 3);

    mses = zeros(numPalettes, 1);
    paletteSize = size(colorPalettes,1)/3;

    % Calculate the MSE for each color palette
    for i = 1:numPalettes


        palette = reshape(colorPalettes(:,i),3,paletteSize)';


        % Calculate the squared Euclidean distance for each pixel to find the closest color
        distances = pdist2(mat, palette, 'euclidean').^2;


        % Find the index of the closest color for each pixel
        [~, closestIndices] = min(distances, [], 2);

        % Generate the new image with colors from the palette
        newImage = palette(closestIndices, :);

        % Reshape the new image to match the original image dimensions
        newImage = reshape(newImage, [], 3);

        % Calculate ate the mean squared error between the original and new image
        mses(i) = mean((double(image(:)) - double(newImage(:))).^2);
    end
end

