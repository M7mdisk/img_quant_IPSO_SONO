function paletteToImg = paletteToImg(colorPalette, img,outPath)

    [height, width,colors] = size(img);

    matrix = double(img);

    if size(matrix, 3) == 4
        matrix = matrix(:, :, 1:3);
    end

    % Normalize the matrix to the range [0, 1]


    mat = reshape(matrix, height * width,3);
    palette = reshape(colorPalette,3,[])';

    % Calculate the squared Euclidean distance for each pixel to find the closest color
    distances = pdist2(mat, palette, 'euclidean').^2;

    % Find the index of the closest color for each pixel
    [~, closestIndices] = min(distances, [], 2);

    % Generate the new image with colors from the palette
    newImage = palette(closestIndices, :);
    % Reshape the new image to match the original image dimensions

    newImage = reshape(newImage, height,width, 3);
    imwrite(newImage/255, outPath);
end

