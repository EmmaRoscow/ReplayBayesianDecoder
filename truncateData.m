
function [start_index, end_index] = truncateData(filename, Timestamps, Xpos, Ypos)

    % Specify start and end of maze exploration, based on stored
    % information (if available) or user visual inspection
    
    switch filename
        
        case 'Achilles_10252013_sessInfo'
            start_index = 4289;
            end_index = length(Timestamps);
            
        case 'Buddy_06272013_sessInfo.mat'
            start_index = 597;
            end_index = length(Timestamps);
            
        otherwise
            % Display data
            f = figure;
            subplot(3, 1, 1); plot(Xpos, Ypos, '.k'); xlabel('Xpos'); ylabel('Ypos')
            title({'Please choose indices where maze exploration starts and ends', '(these should be integers). Then press any key.'})
            subplot(3, 1, 2); plot(Xpos, '.k'); ylabel('Xpos'); title(['Length of data: ' num2str(length(Timestamps))])
            subplot(3, 1, 3); plot(Ypos, '.k'); ylabel('Ypos')
            pause
            
            % Ask user for start and end index
            idx = inputdlg({'Please give index where maze exploration starts', 'Please give index where maze exploration ends'});
            evalc('start_index = str2num(idx{1});');
            evalc('end_index = str2num(idx{2});');
            close(f)
            
            % Confirm for user what this truncation looks like
            figure
            
            subplot(3, 1, 1); hold on
            plot(Xpos(start_index:end_index), Ypos(start_index:end_index), '.k');
            plot(Xpos(1:start_index-1), Ypos(1:start_index-1), '.b');
            plot(Xpos(end_index+1:end), Ypos(end_index+1:end), '.b');
            xlabel('Xpos'); ylabel('Ypos')
            
            subplot(3, 1, 2); hold on
            plot(start_index:end_index, Xpos(start_index:end_index), '.k');
            plot(1:start_index-1, Xpos(1:start_index-1), '.b');
            plot(end_index+1:length(Xpos), Xpos(end_index+1:end), '.b');
            ylabel('Xpos');
            
            subplot(3, 1, 3); hold on
            plot(start_index:end_index, Ypos(start_index:end_index), '.k');
            plot(1:start_index-1, Ypos(1:start_index-1), '.b');
            plot(end_index+1:length(Ypos), Ypos(end_index+1:end), '.b');
            ylabel('Ypos')
            
    end
    
end