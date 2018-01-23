classdef BehaviorData < hgsetget %& dynamicprops
    
    properties
        basic_info
        trial_dataset
        shuffle_dataset
    end
    
    properties
        shuffle_dataset_count = 10;
    end
    
    events
    end
    
    methods
        %% constructor
        function obj = BehaviorData(trial_file, lick_file)
            if ~nargin
                fprintf('Select TRIAL file....\n');
                [fname,fdir] = uigetfile('*.*','MultiSelect','off');
                cd(fdir);
                trial_filename = [fdir,fname];
                fprintf('Select LICK file....\n');
                [fname,fdir] = uigetfile('*.*','MultiSelect','off');
                cd(fdir);
                lick_filename = [fdir,fname];
            end
            
            trial_data = load(trial_filename);
            lick_data = load(lick_filename);
            
            obj.basic_info.trial_filename = trial_filename;
            obj.basic_info.lick_filename = lick_filename;
            
            obj.basic_info.sound_duration = unique(trial_data(:,5))/1000; % convert to second
            if length(obj.basic_info.sound_duration)>1
                fprintf('Inconsistent sound duration\n');
            end
            
            obj.basic_info.reaction_duration = unique(trial_data(:,6))/1000; % convert to second
            if length(obj.basic_info.reaction_duration)>1
                fprintf('Inconsistent rection duration\n');
            end
            
            obj.basic_info.total_trial_count = size(trial_data,1);
            obj.basic_info.total_lick_count = size(lick_data,1);
            
            obj.basic_info.sound_type = unique(trial_data(:,3));
            
            obj.basic_info.trial_time = trial_data(:,1);
            obj.basic_info.trial_sound = trial_data(:,3);
            obj.basic_info.trial_ITI = trial_data(:,7)/1000+1.5; % convert to second and add the "pre-trial"
            obj.basic_info.trial_length = obj.basic_info.sound_duration+obj.basic_info.reaction_duration+obj.basic_info.trial_ITI;
            
            obj.basic_info.lick_time = lick_data(:,1);
            
            obj.basic_info.trace_time = [];
            
            for i=1:obj.basic_info.total_trial_count
                obj.trial_dataset.trial(i).sound = obj.basic_info.trial_sound(i);
                obj.trial_dataset.trial(i).length = obj.basic_info.trial_length(i);
                obj.trial_dataset.trial(i).start_time = obj.basic_info.trial_time(i);
                obj.trial_dataset.trial(i).lick_time = obj.basic_info.lick_time-obj.trial_dataset.trial(i).start_time;
                
                if find(obj.trial_dataset.trial(i).lick_time>0 & obj.trial_dataset.trial(i).lick_time<(obj.basic_info.sound_duration+obj.basic_info.reaction_duration),1)
                    obj.trial_dataset.trial(i).water_delivered = 1;
                else
                    obj.trial_dataset.trial(i).water_delivered = 0;
                end
                
                if find(obj.trial_dataset.trial(i).lick_time>0 & obj.trial_dataset.trial(i).lick_time<obj.basic_info.sound_duration,1)
                    obj.trial_dataset.trial(i).lick_during_sound = 1;
                else
                    obj.trial_dataset.trial(i).lick_during_sound = 0;
                end
                
                obj.trial_dataset.trial(i).first_lick_idx = find(obj.trial_dataset.trial(i).lick_time>0 & obj.trial_dataset.trial(i).lick_time<obj.trial_dataset.trial(i).length,1);
                obj.trial_dataset.trial(i).first_lick_time = obj.trial_dataset.trial(i).lick_time(obj.trial_dataset.trial(i).first_lick_idx);
                obj.trial_dataset.trial(i).lick_with_water_idx = find(obj.trial_dataset.trial(i).lick_time>0 & obj.trial_dataset.trial(i).lick_time<(obj.basic_info.sound_duration+obj.basic_info.reaction_duration),1);
                obj.trial_dataset.trial(i).lick_with_water_time = obj.trial_dataset.trial(i).lick_time(obj.trial_dataset.trial(i).lick_with_water_idx);
                %ITI lick is defined as the last lick occured after 1 second of reaction window and before 1 second of the enxt trial 
                obj.trial_dataset.trial(i).ITI_lick_idx = find(obj.trial_dataset.trial(i).lick_time>(obj.basic_info.sound_duration+obj.basic_info.reaction_duration+1) & obj.trial_dataset.trial(i).lick_time<(obj.trial_dataset.trial(i).length-1),1,'last');
                obj.trial_dataset.trial(i).ITI_lick_time = obj.trial_dataset.trial(i).lick_time(obj.trial_dataset.trial(i).ITI_lick_idx);
                
                obj.trial_dataset.trial(i).trace_time = [];
            end
            
            obj.trial_dataset.correct_rate = mean(cat(1,obj.trial_dataset.trial.lick_during_sound)); % 'correct' means lick during the sound
            
%             fprintf(['Generating ',num2str(obj.shuffle_dataset_count),' shuffle datasets....\n']);
            for r=1:obj.shuffle_dataset_count
                if mod(r,10)==0
                    fprintf(['Generating shuffle datasets ',num2str(r),'/',num2str(obj.shuffle_dataset_count),'\n']);
                end
                
                shuffle_idx = randperm(obj.basic_info.total_trial_count);
                shuffle_trial_length = obj.basic_info.trial_length(shuffle_idx);
                shuffle_trial_sound = obj.basic_info.trial_sound(shuffle_idx);
                
                for i=1:obj.basic_info.total_trial_count
                    obj.shuffle_dataset(r).trial(i).sound = shuffle_trial_sound(i);
                    obj.shuffle_dataset(r).trial(i).length = shuffle_trial_length(i);
                    
                    if i==1
                        obj.shuffle_dataset(r).trial(i).start_time = obj.basic_info.trial_time(1);
                    else
                        obj.shuffle_dataset(r).trial(i).start_time = obj.shuffle_dataset(r).trial(i-1).start_time+obj.shuffle_dataset(r).trial(i-1).length;
                    end
                    
                    obj.shuffle_dataset(r).trial(i).lick_time = obj.basic_info.lick_time-obj.shuffle_dataset(r).trial(i).start_time;
                    
                    if find(obj.shuffle_dataset(r).trial(i).lick_time>0 & obj.shuffle_dataset(r).trial(i).lick_time<(obj.basic_info.sound_duration+obj.basic_info.reaction_duration))
                        obj.shuffle_dataset(r).trial(i).water_delivered = 1;
                    else
                        obj.shuffle_dataset(r).trial(i).water_delivered = 0;
                    end

                    if find(obj.shuffle_dataset(r).trial(i).lick_time>0 & obj.shuffle_dataset(r).trial(i).lick_time<obj.basic_info.sound_duration)
                        obj.shuffle_dataset(r).trial(i).lick_during_sound = 1;
                    else
                        obj.shuffle_dataset(r).trial(i).lick_during_sound = 0;
                    end
                    
                    obj.shuffle_dataset(r).correct_rate = mean(cat(1,obj.shuffle_dataset(r).trial.lick_during_sound)); % 'correct' means lick during the sound
                    
                    obj.shuffle_dataset(r).trial(i).trace_time = [];
                end
                
            end
            
            
            
            
        end
        
        %% plot
        
        function plot_lick(obj, data_type, window_range)
            if ~ismember(data_type, {'exp-all', 'exp-sound', 'shuffle'})
                fprintf('data_type has to be all or sound\n');
            else
                switch data_type
                    case 'exp-all'
                        fprintf('Plot experiment dataset, combining both sounds.\n');
                        figure
                        hold on;                   
                        for i=1:obj.basic_info.total_trial_count
                            data = obj.trial_dataset.trial(i).lick_time(obj.trial_dataset.trial(i).lick_time>=window_range(1) & obj.trial_dataset.trial(i).lick_time<=window_range(2));
                            if ~isempty(data)
                                plot(data,i,'b.')
                            end
                        end
                        hold off;
                        xlim(window_range);
                        ylim([0 obj.basic_info.total_trial_count+1]);
                        set(gca,'YDir','reverse');
                    case 'exp-sound'
                        fprintf(['Plot experiment dataset, seperate ', num2str(size(obj.basic_info.sound_type,1)), ' sounds.\n']);
                        for s=1:size(obj.basic_info.sound_type,1)
                            sound_idx = find([obj.trial_dataset.trial.sound]==obj.basic_info.sound_type(s));
                            figure;
                            hold on;
                            for i=1:length(sound_idx)
                                data = obj.trial_dataset.trial(sound_idx(i)).lick_time(obj.trial_dataset.trial(sound_idx(i)).lick_time>=window_range(1) & obj.trial_dataset.trial(sound_idx(i)).lick_time<=window_range(2));
                                if ~isempty(data)
                                    plot(data,i,'b.')
                                end
                            end
                            hold off;
                            title(['Sound ',num2str(obj.basic_info.sound_type(s))]);
                            xlim(window_range);
                            ylim([0 length(sound_idx)+1]);
                            set(gca,'YDir','reverse');
                        end
                    case 'shuffle'
                        fprintf('Plot single shuffle dataset.\n');
                        shuffle_dataset_id = input(['Select the shuffle dataset (1-',num2str(obj.shuffle_dataset_count),'):']);
                        figure
                        hold on;                   
                        for i=1:obj.basic_info.total_trial_count
                            data = obj.shuffle_dataset(shuffle_dataset_id).trial(i).lick_time(obj.shuffle_dataset(shuffle_dataset_id).trial(i).lick_time>=window_range(1) & obj.shuffle_dataset(shuffle_dataset_id).trial(i).lick_time<=window_range(2));
                            if ~isempty(data)
                                plot(data,i,'b.')
                            end
                        end
                        hold off;
                        title(['Shuffle dataset ',num2str(shuffle_dataset_id)]);
                        xlim(window_range);
                        ylim([0 obj.basic_info.total_trial_count+1]);
                        set(gca,'YDir','reverse');
                end
            end
        end
        
        function plot_lick_frequency(obj, data_type, window_range, bin_size)
            if ~ismember(data_type, {'exp-all', 'exp-sound', 'shuffle'})
                fprintf('data_type has to be all or sound\n');
            else
                switch data_type
                    case 'exp-all'
                        fprintf('Plot experiment dataset, combining both sounds.\n');
                        data = cat(1,obj.trial_dataset.trial.lick_time);
                        data = data(all_lick_time>=window_range(1) & all_lick_time<=window_range(2));
                        hist_x = window_range(1):bin_size:window_range(2);
                        bincounts = histc(data,hist_x);
                        bincounts = bincounts/obj.basic_info.total_trial_count/bin_size;
                        figure;
                        bar(hist_x,bincounts,'histc');
                        xlim(window_range);
                        ylabel('lick per second');
                    case 'exp-sound'
                        fprintf(['Plot experiment dataset, seperate ', num2str(size(obj.basic_info.sound_type,1)), ' sounds.\n']);

                        for s=1:size(obj.basic_info.sound_type,1)
                            data = cat(1,obj.trial_dataset.trial([obj.trial_dataset.trial.sound]==obj.basic_info.sound_type(s)).lick_time);
                            data = data(all_lick_time>=window_range(1) & all_lick_time<=window_range(2));
                            hist_x = window_range(1):bin_size:window_range(2);
                            bincounts = histc(data,hist_x);
                            bincounts = bincounts/obj.basic_info.total_trial_count/bin_size;
                            figure;
                            bar(hist_x,bincounts,'histc');
                            title(['Sound ',num2str(obj.basic_info.sound_type(s))]);
                            xlim(window_range);
                            ylabel('lick per second');
                        end
                    case 'shuffle'
                        fprintf('Plot single shuffle dataset.\n');
                        shuffle_dataset_id = input(['Select the shuffle dataset (1-',num2str(obj.shuffle_dataset_count),'):']);
                        data = cat(1,obj.shuffle_dataset(shuffle_dataset_id).trial.lick_time);
                        hist_x = window_range(1):bin_size:window_range(2);
                        bincounts = histc(data,hist_x);
                        bincounts = bincounts/obj.basic_info.total_trial_count/bin_size;
                        figure;
                        bar(hist_x,bincounts,'histc');
                        title(['Shuffle dataset ',num2str(shuffle_dataset_id)]);
                        xlim(window_range);
                        ylabel('lick per second');      
                end
            end
                
        end
        
        
        function plot_correct_rate_over_time(obj, data_type, window_size)
            if ~ismember(data_type, {'exp-all', 'exp-sound', 'exp-sound', 'shuffle', 'shuffle-all-dataset', 'both', 'both-all-dataset'})
                fprintf('data_type has to be exp, random, or both\n');
            else
                switch data_type
                    case 'exp-all'
                        fprintf('Plot experiment dataset, combining both sounds.\n');
                        data = cat(1,obj.trial_dataset.trial.lick_during_sound);
                        out = average_over_moving_window(data, window_size);
                        plot(out);
                    case 'exp-sound'
                        fprintf(['Plot experiment dataset, seperate ', num2str(size(obj.basic_info.sound_type,1)), ' sounds.\n']);

                        for s=1:size(obj.basic_info.sound_type,1)
                            data = cat(1,obj.trial_dataset.trial([obj.trial_dataset.trial.sound]==obj.basic_info.sound_type(s)).lick_during_sound);
                            out = average_over_moving_window(data, window_size);
                            figure;
                            plot(out);
                            title(['Sound ',num2str(obj.basic_info.sound_type(s))]);
                            ylim([0 1.05]);
                        end
 
                    case 'shuffle'
                        fprintf('Plot single shuffle dataset.\n');
                        shuffle_dataset_id = input(['Select the shuffle dataset (1-',num2str(obj.shuffle_dataset_count),'):']);
                        data = cat(1,obj.shuffle_dataset(shuffle_dataset_id).trial.lick_during_sound);
                        out = average_over_moving_window(data, window_size);
                        figure;
                        plot(out);
                        title(['Shuffle dataset ',num2str(shuffle_dataset_id)]);
                        ylim([0 1.05]);
                        
                    case 'shuffle-all-dataset'
                        fprintf('Plot all shuffle dataset.\n');
                        out = zeros(obj.basic_info.total_trial_count-window_size+1,obj.shuffle_dataset_count);
                        for i=1:obj.shuffle_dataset_count
                            data = cat(1,obj.shuffle_dataset(i).trial.lick_during_sound);
                            out(:,i) = average_over_moving_window(data, window_size);
                        end
                        plot(out,'b.');
                    case 'both'
                    case 'both-all-dataset'
                        data = cat(1,obj.trial_dataset.trial.lick_during_sound);
                        out = average_over_moving_window(data, window_size);
                        plot(out,'-r.');
                        hold on;
                        out = zeros(obj.basic_info.total_trial_count-window_size+1,obj.shuffle_dataset_count);
                        for i=1:obj.shuffle_dataset_count
                            data = cat(1,obj.shuffle_dataset(i).trial.lick_during_sound);
                            out(:,i) = average_over_moving_window(data, window_size);
                        end
                        plot(out,'b.');
                        hold off;
                end

                ylim([0 1.05]);
            end
        end
        
        function plot_random_correct_rate_distribution(obj)
            hist_x = 0:0.05:1;
            bincounts = histc(cat(1,obj.shuffle_dataset.correct_rate),hist_x);
            bar(hist_x,bincounts,'histc');
            xlim([0 1]);
            
            line([obj.trial_dataset.correct_rate obj.trial_dataset.correct_rate],[0 max(bincounts)],'color','r')
            
        end
        
    end
end

function [out] = average_over_moving_window(data, window_size)
    if window_size>length(data)
        warning('Incorrect moving-window size');
    else
        out = zeros(length(data)-window_size+1,1);
        for i=1:length(out)
            out(i) = mean(data(i:i+window_size-1));
        end
    end
        
end


