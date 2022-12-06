function output_imgs = new_LoadImageStacks( inputPath, dims, sub_dir, useGPU )
%LOADIMAGESTACKS Load image stacks for each round

    % Suppress all warnings 
    warning('off','all');
    
    % Set output_imgs
    output_imgs = zeros(dims, 'uint8');
    
    % Get directories containing all images 
    dirs = dir(strcat(inputPath, 'round*'));
    Nround = numel(dirs);
    Nchannel = 4;
    
    parfor r=1:Nround
%     for r=1:Nround
        tic
        fprintf('Loading round %d...', r);
        
        curr_dir = dirs(r).name;
        
        curr_files = dir(fullfile(inputPath, curr_dir, sub_dir, '*.tif'));

        %Nchannel = numel(curr_files);
       
        
        % Load all channels
        for c=1:Nchannel 
            curr_path = strcat(curr_files(c).folder, '/', curr_files(c).name);

            curr_img = new_LoadMultipageTiff(curr_path, 'unit8', 'uint8', useGPU);
            
%             if isempty(zrange)
%                 round_img(:,:,:,c) = curr_img;
%             else
%                 round_img(:,:,:,c) = curr_img(:,:,zrange);
%             end
            
            output_imgs(:,:,:,c,r) = curr_img;
        end
        fprintf(sprintf('[time = %.2f s]\n', toc));
    end

    % Collapse to common sized array
%     maxX = 1E10; maxY = 1E10;
%     for r=1:Nround
%         curr_round = cell_imgs{r};
%        [currX,currY,~] = size(curr_round); 
%        if currX < maxX
%            maxX = currX;
%        end
%        if currY < maxY
%            maxY = currY;
%        end
%     end
    
%     dims = [maxX maxY size(curr_round, 3) size(curr_round, 4) Nround];
% 
%     if useGPU
%         output_imgs = zeros(dims, 'uint8', 'gpuArray');
%     else
%         output_imgs = zeros(dims, 'uint8');
%     end
    
    % Show message for re-sizing
%     fprintf('Collapsed to size %d by %d by %d\n', maxX, maxY, size(curr_round, 3));
    
%     for r=1:Nround
%         fprintf('Collapsing round %d\n', r);
%         curr_round = cell_imgs{r};
%         output_imgs(:,:,:,:,r) = curr_round(1:maxX, 1:maxY, :, :);
%     end

end

