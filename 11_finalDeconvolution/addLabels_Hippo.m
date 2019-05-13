function labels = addLabels_Hippo(coordinates)
% Create ZBrainAtlas label matrix for coordinates in ZBrainAtlas space

    % display
    disp('getting labels');

    maskdbPath = '/home/ljp/Science/Hippolyte/MaskDatabase.mat'; % TODO insert in focus

    % Load ZBrainAtlas (from https://engertlab.fas.harvard.edu/Z-Brain/#/downloads)
    load(maskdbPath, 'DateCreated', 'MaskDatabase', 'MaskDatabaseNames', 'MaskDatabaseOutlines', 'height', 'width', 'Zs')

    % This yields:
    % Dimensions: height × width × Zs
    % Region names: MaskDatabaseNames
    % Sparse representation of regions: MaskDatabase
    % Sparse representation of region outlines: MaskDatabaseOutlines
    % Documentation: resolution: x/y/z = 0.798/0.798/2um

    % Resolution from ZBrainAtlas doc in um
    resolution = [0.798, 0.798, 2];

    % Find closest voxel per coordinate
    % Assuming(!) that coordinates is already in the same space as the
    % ZBrainAtlas, one can find the closest voxel center per coordinate (i.e.
    % the voxel that contains the coordinate).
    grid_coordinates = round(coordinates * 1000 ./ resolution);
    % Calculates the data coordinates expressed in the ZBrainAtlas grid (go to um and apply resolution transformation)

    % Flip coordinates (determined after visual inspection 11/06/2018)
    grid_coordinates = [grid_coordinates(:, 2), grid_coordinates(:,1), grid_coordinates(:, 3)]; % reorientate x and y axis
    grid_coordinates(:, 1) = height - grid_coordinates(:, 1); % swap new x coordinates left to right

    % optionally plot brains
    % plotBrains()

    % Labeling
    labels = zeros(length(coordinates), length(MaskDatabaseNames)); % must be full matrix instead of sparse matrix, because of transfer to python
    count_out_of_bound = 0;

    for iNeuron = 1:length(coordinates) % can this be (more efficient) without a loop?
        % Check if within atlas-space, fetch label
        if (grid_coordinates(iNeuron, 3) <= Zs) && (grid_coordinates(iNeuron, 2) <= width) && (grid_coordinates(iNeuron, 1) <= height)
            labels(iNeuron, :) = MaskDatabase(sub2ind([height width Zs], grid_coordinates(iNeuron,1), grid_coordinates(iNeuron,2), grid_coordinates(iNeuron,3)),:); % sparse array is put into full matrix
        else
            count_out_of_bound = count_out_of_bound + 1; 
        end
    end

    fprintf('out of bounds : %d\n', count_out_of_bound);

    % Export labels
    % labels now contains the labels of all neurons, add this to the hdf5 file
    % in order to Fishualize.

    function plotBrains()
        % Plotting (default should be 0)
        % enable to manually inspect the outline plot of ZBrainAtlas with
        % grid_coordinates in a three 2D scatter plots

        % Create ZBrainAtlas coordinates for plotting (NB: Computationally heavy!)
        % one needs the findND.m function to go with this, link:
        % https://uk.mathworks.com/matlabcentral/fileexchange/64383-findnd
        % I used v1.0
        dummy_coordinates = ones(height, width, Zs); %h/w/Z as defined by ZBrain
        [atl_x, atl_y, atl_z] = findND(dummy_coordinates); % Use custom 3D find function to catch coordinates
        clear dummy_coordinates; % clear memory
        % resolution = [0.798, 0.798, 2];
        resolution = [1, 1, 1];
        atlas_coordinates =  [(resolution(1) *atl_x)  (resolution(2) *atl_y) (resolution(3) * atl_z)]; % resize to resolutions (Z Brain Atlas Github): x/y/z = 0.798/0.798/2um.
        clear atl_x; clear atl_y; clear atl_z;  %clear memory
        % atlas_coordinates is now defined in um
        max_atlas_coords = max(atlas_coordinates);
        min_atlas_coords = min(atlas_coordinates);

        % Plot atlas outline

        close('all')
        skip = 100;
        test = sum(MaskDatabaseOutlines');
        testfind = find(test);
        atlas_fish = atlas_coordinates(testfind, :);
        figure(4)
        scatter(grid_coordinates(:,1), grid_coordinates(:,2),'b')
        hold('on')
        scatter(atlas_fish(1:skip:end,1), atlas_fish(1:skip:end,2),'r')
        xlabel('x'); ylabel('y'); grid on;

        figure(5)
        scatter(grid_coordinates(:,2), grid_coordinates(:,3),'b')
        hold('on')
        scatter(atlas_fish(1:skip:end,2), atlas_fish(1:skip:end,3),'r')
        xlabel('y'); ylabel('z'); grid on;

        figure(6)
        scatter(grid_coordinates(:,1), grid_coordinates(:,3),'b')
        hold('on')
        scatter(atlas_fish(1:skip:end,1), atlas_fish(1:skip:end,3),'r')
        xlabel('x'); ylabel('z'); grid on;

    end

end
