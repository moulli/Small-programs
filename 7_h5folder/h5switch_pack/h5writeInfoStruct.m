function information = h5writeInfoStruct(F, paramsIN)

%% Function that interactively creates the three structures for h5switch.
%
%  Some information is given be the user, and then the function fills the
%  blanks in the three structures by asking directly the user in the
%  console.
%
%
%% Parameters:
%  
%  --F: the focus of interest.
%  --paramsIN: if there are parameters that are already given, there are
%    here.
%
%
%% Outputs
%
%  --information: structure for h5switch.


    
    %% Getting date and run:
    
    filetemp = F.date;
    runtemp = F.run;
    
    
    
    %% Creating information structure:
    
    information = struct;
    information.date = filetemp;
    information.run = str2double(runtemp(end-2:end));
    
    
    
    %% Fields completed:
    
    ncomplete = zeros(10, 1);
    while sum(ncomplete) ~= 10
    
    
    
    %% Description:
    
        while ncomplete(1) == 0
            fprintf('\n');
            if isfield(paramsIN, 'description')
                information.description = paramsIN.description;
                ncomplete(1) = 1;
            elseif isfield(paramsIN, 'autodescription')
                information.autodescription = paramsIN.autodescription;
                ncomplete(1) = 1;
            else
                prompt = 'Please provide a description. [N] to skip, [Auto] for autodescritpion, [Exit] to stop process: ';
                instr = input(prompt, 's');
                switch lower(instr)
                    case 'exit'
                        error('Program aborted.')
                    case 'n'
                        information.autodescription = false;
                        ncomplete(1) = 1;
                    case 'auto'
                        information.autodescription= true;
                        ncomplete(1) = 1;
                    otherwise
                        if isempty(instr)
                            continue
                        end
                        information.autodescription = false;
                        information.description = instr;
                        ncomplete(1) = 1;
                end
            end
        end



        %% Acquisition information:

        % Rate:
        while ncomplete(2) == 0
            fprintf('\n');
            if isfield(paramsIN, 'rate')
                information.rate = paramsIN.rate;
                ncomplete(2) = 1;
            else
                prompt = 'Please provide acquisition rate (Hz). [Exit] to stop process: ';
                instr = input(prompt, 's');
                switch lower(instr)
                    case 'exit'
                        error('Program aborted.')
                    otherwise
                        instr = str2double(instr);
                        if ~isnan(instr)
                            information.rate = instr;
                            ncomplete(2) = 1;
                        else
                            fprintf('Acquisition rate must be a number. \n');
                        end
                end   
            end
        end

        % Layers:
        while ncomplete(3) == 0
            if isfield(paramsIN, 'layers')
                information.layers = paramsIN.layers;
                ncomplete(3) = 1;
            else
                prompt = 'Please provide number of layers. [Exit] to stop process: ';
                instr = input(prompt, 's');
                switch lower(instr)
                    case 'exit'
                        error('Program aborted.')
                    otherwise
                        instr = str2double(instr);
                        if ~isnan(instr) && (instr == floor(instr))
                            information.layers = instr;
                            ncomplete(3) = 1;
                        else
                            fprintf('Number of layers must be an integer. \n');
                        end
                end   
            end
        end

        % Increment:
        while ncomplete(4) == 0
            if isfield(paramsIN, 'increment')
                information.increment = paramsIN.increment;
                ncomplete(4) = 1;
            else
                prompt = 'Please provide interlayer distance (µm). [Exit] to stop process: ';
                instr = input(prompt, 's');
                switch lower(instr)
                    case 'exit'
                        error('Program aborted.')
                    otherwise
                        instr = str2double(instr);
                        if ~isnan(instr)
                            information.increment = instr;
                            ncomplete(4) = 1;
                        else
                            fprintf('Interlayer distance must be a number. \n');
                        end
                end  
            end
        end



        %% Fish information:

        % Fish line:
        while ncomplete(5) == 0
            fprintf('\n');
            if isfield(paramsIN, 'line')
                information.line = paramsIN.line;
                ncomplete(5) = 1;
            else
                prompt = 'Please provide fish line. [Exit] to stop process: ';
                instr = input(prompt, 's');
                switch lower(instr)
                    case 'exit'
                        error('Program aborted.')
                    otherwise
                        if isempty(instr)
                            continue
                        end
                        information.line = instr;
                        ncomplete(5) = 1;
                end   
            end
        end

        % Fish age:
        while ncomplete(6) == 0
            if isfield(paramsIN, 'age')
                information.age = paramsIN.age;
                ncomplete(6) = 1;
            else
                prompt = 'Please provide fish age (dpf). [Exit] to stop process: ';
                instr = input(prompt, 's');
                switch lower(instr)
                    case 'exit'
                        error('Program aborted.')
                    otherwise
                        instr = str2double(instr);
                        if ~isnan(instr) && (instr == floor(instr))
                            information.age = instr;
                            ncomplete(6) = 1;
                        else
                            fprintf('Fish age must be an integer. \n');
                        end
                end  
            end
        end

        % Fish ID:
        while ncomplete(7) == 0
            if isfield(paramsIN, 'ID')
                information.ID = paramsIN.ID;
                ncomplete(7) = 1;
            elseif isfield(paramsIN, 'autoID')
                information.autoID = paramsIN.autoID;
                ncomplete(7) = 1;
            else
                prompt = 'Please provide the fish ID. [Auto] for auto ID, [Exit] to stop process: ';
                instr = input(prompt, 's');
                switch lower(instr)
                    case 'exit'
                        error('Program aborted.')
                    case 'auto'
                        information.autoID= true;
                        ncomplete(7) = 1;
                    otherwise
                        if isempty(instr)
                            continue
                        end
                        information.autoID = false;
                        information.ID = instr;
                        ncomplete(7) = 1;
                end
            end
        end



        %% Stimuli:

        stim = struct;
        while ncomplete(8) == 0
            fprintf('\n');
%             prompt = 'Please provide number of stimuli (can be 0). [Exit] to stop process: ';
%             instr = input(prompt, 's');
%             instr = str2double(instr);
            instr = 1;
            if isnan(instr) || (instr ~= floor(instr))
                fprintf('Number of stimuli must be an integer. \n');
            else
                stim.name = cell(1, instr);
                stim.sensorytype = cell(1, instr);
                stim.stimtype = cell(1, instr);
                stim.frequency = cell(1, instr);
                for k = 1:instr
                    fprintf('Please provide information for stimulus number %.0f.\n', k);
                    while true
%                         prompt2 = 'Please provide name of stimulus:';
%                         instr2 = input(prompt2, 's');
                        instr2 = 'vestibular1';
                        if ~isempty(instr2)
                            stim.name{k} = instr2;
                            break
                        end
                    end
                    while true
%                         prompt2 = 'Please provide sensory type (e.g. [vestibular], [thermotaxis]...):';
%                         instr2 = input(prompt2, 's');
                        instr2 = 'vestibular';
                        if ~isempty(instr2)
                            stim.sensorytype{k} = instr2;
                            break
                        end
                    end
                    while true
                        prompt2 = 'Please provide stimulus type (e.g. [sine], [step]...):';
                        instr2 = input(prompt2, 's');
                        if ~isempty(instr2)
                            stim.stimtype{k} = instr2;
                            break
                        end
                    end
                    while true
                        prompt2 = 'Please provide frequency (Hz). [0] if aperiodic:';
                        instr2 = input(prompt2, 's');
                        instr2 = str2double(instr2);
                        if isnan(instr2)
                            fprintf('Frequency must be a number. \n');
                            continue
                        else
                            stim.frequency{k} = instr2;
                        end
                        break
                    end
                end
                ncomplete(8) = 1;
            end  
        end
        information.stimulus = stim;



        %% Behaviours:

        behav = struct;
        while ncomplete(9) == 0
            fprintf('\n');
%             prompt = 'Please provide number of behaviours (can be 0). [Exit] to stop process: ';
%             instr = input(prompt, 's');
%             instr = str2double(instr);
            instr = 0;
            if isnan(instr) || (instr ~= floor(instr))
                fprintf('Please provide number of behaviours as an integer. \n');
            else
                behav.name = cell(1, instr);
                behav.type = cell(1, instr);
                behav.frequency = cell(1, instr);
                for k = 1:instr
                    fprintf('Please provide information for behaviour number %.0f.\n', k);
                    while true
                        prompt2 = 'Please provide name of behaviour:';
                        instr2 = input(prompt2, 's');
                        if ~isempty(instr2)
                            behav.name{k} = instr2;
                            break
                        end
                    end
                    while true
                        prompt2 = 'Please provide type (e.g. [eye tracking], [tail tracking]...):';
                        instr2 = input(prompt2, 's');
                        if ~isempty(instr2)
                            behav.type{k} = instr2;
                            break
                        end
                    end
                    while true
                        prompt2 = 'Please provide frequency (Hz). [0] if aperiodic:';
                        instr2 = input(prompt2, 's');
                        instr2 = str2double(instr2);
                        if isnan(instr2)
                            fprintf('Frequency must be a number. \n');
                            continue
                        else
                            behav.frequency{k} = instr2;
                        end
                        break
                    end
                end
                ncomplete(9) = 1;
            end  
        end
        information.behaviour = behav;



        %% Checking if information is correct:

        information
        prompt = 'Are these fields correct? [Y] or [N]:';
        instr = input(prompt, 's');
        switch lower(instr)
            case 'y'
                break
            case 'n'
                while true
                    prompt2 = 'Which field do you want to change? [None] to exit:';
                    instr2 = input(prompt2, 's');
                    switch lower(instr2)
                        case 'none'
                            break
                        otherwise
                            if ~isfield(information, instr2)
                                fprintf('Please enter an actual field. \n');
                            else
                                fields = {'description', 'rate', 'layers', 'increment', 'line', 'age', 'ID', 'stimulus', 'behaviour'};
                                for z = 1:length(fields)
                                    if isequal(instr2, fields{z})
                                        ncomplete(z) = 0;
                                        break
                                    end
                                end
                                break
                            end
                    end
                end
                continue
            otherwise
                fprintf('Please enter [Y] or [N]. \n');
                continue
        end
        
    end
    
    
end