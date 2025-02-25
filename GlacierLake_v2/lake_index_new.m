function [bot_lake_ind, bot_lake_depth, top_lake_ind, top_lake_depth, slush_lid_top, slush_lid_bot] =...
    lake_index_new(lambda, ice_lim, total_grid_num, grid_profile, tt, t_step, slush_lid_threshold, slush_lid_num)
%GlacierLake function to index lakes. This is the only subfunction that
%hasn't had a major spring clean for the github upload due to length,
%complexity and my lack of time

    %ice_lim = 2; %number of cells that need to not be ice in order for ice lid to be treated as seperate
    slush_in_lid_max = 3; %maximum number of slush cells allowed inside lid comprised of ice
    slush_track = 0; %tracker

    %this code will handle 2 ice lids and two lakes. If only one lake is found
    %it is referred to with the 'bot' prefix. If two lakes are found the lowest
    %is referred to with the 'bot' prefix and the top with the 'top' prefix
    
    %set empty variables to avoid errors and allow correct paths to be taken
    top_ice_lid_top = [];
    top_ice_lid_bot = [];
    top_lake_ind = [];
    bot_ice_lid_top = [];
    bot_ice_lid_bot = [];
    bot_lake_ind = [];
    ice_substrate_top = [];
    slush_lid_top = [];
    slush_lid_bot = [];
    bot_lake_bot = [];

    %calculate top lake depth and location. Using first slush cell rather than first
    %ice cell
    water_ind = lambda >= 1;
    slush_ind = lambda > 0 & lambda < 1;
    ice_ind = lambda == 0;

    %this gives the last ice cell upwards from the base
    ice_substrate_top = total_grid_num - (min(strfind(fliplr(ice_ind'), [1 0])) - 1);

    if isempty(ice_substrate_top)
        %warning('No lake found.'); fprintf('%.2f', (tt*t_step)/24)
        ice_substrate_top = 1;
    end

    %if there is ice above the substrate then there is at least one ice lid to
    %account for
    if any(ice_ind(1:ice_substrate_top - 1) == 1)

        %this finds the bottom ice cell of the lowest ice lid. Looks complicated
        %but is robust
        ice_ind_sub = find(ice_ind);
        ice_ind_sub_ud = flipud(ice_ind_sub);
        bot_ice_lid_bot = ice_ind_sub_ud(find(diff(ice_ind_sub_ud) < - ice_lim, 1, 'first') + 1);%盖子底部位置

        %this finds the top of the lowest ice lid
        bot_ice_lid_top = bot_ice_lid_bot - strfind(flipud(ice_ind(1:bot_ice_lid_bot))', [1 0]) + 1;

        %if the above method gives an empty value then the bottom ice lid is the
        %only ice lid and therefore also the top ice lid so another method needs to
        %be used
        if isempty(bot_ice_lid_top)
            bot_ice_lid_top = bot_ice_lid_bot - find(flipud(ice_ind(1:bot_ice_lid_bot)) == 1, 1, 'last') + 1;
        end

        bot_ice_lid_top = bot_ice_lid_top(1);
    end
        
    %if there is no ice above the top of the ice substrate then the lake index
    %can be obtanied. If ice_substrate_top = 1 then it is the upmost cell and
    %there is therefore no lake
    if ice_substrate_top ~= 1

        %find the lowest water cell of the upmost lake
        bot_lake_bot = total_grid_num - min(strfind(fliplr(water_ind'), [0 1]));

        %this means that surface slush is not accounted for as part of the lake
        %unless there is only slush in which case only the top slush cell is
        %taken as part of the lake
        bot_lake_top = find(water_ind == 1, 1, 'first');

        %if there is no water cell locatable then lake_top and lake_ind = 1
        if isempty(bot_lake_top)
            bot_lake_top = 1;
            bot_lake_ind = 1;

        %if there is water then lake_ind can be calculated
        else    
            bot_lake_ind = (bot_lake_top:bot_lake_bot)';
        end

    %if non of these things it must be an error
    elseif ice_substrate_top == 1
        %no lake
    else
        error('An unacounted for scenario has arisen.'); fprintf('%.2f', (tt*t_step)/24)
    end

    %IF A LID IS PRESENT

    %if there is nothing above the bottom lid
    if bot_ice_lid_top == 1

        %find the lowest water cell of the  lake
        bot_lake_bot = total_grid_num - min(strfind(fliplr(water_ind'), [0 1]));

        %find the top water cell of the lake
        bot_lake_top = find(water_ind == 1, 1, 'first');

        bot_lake_ind = (bot_lake_top:bot_lake_bot)';


    %if there are more cells above bot_ice_lid_top then there is certainly
    %water or slush and possibly ice above the bottom lid that must be
    %accounted for

    %if there is ice and water above the lowest lid
    elseif any(ice_ind(1:bot_ice_lid_top - 1) == 1) &&...
            any(water_ind(1:bot_ice_lid_top) == 1) 

        %this line will find the bottom of the top ice lid
        top_ice_lid_bot = min(ice_ind_sub_ud(find(diff(ice_ind_sub_ud) < - ice_lim, 2, 'first') + 1));

        %this finds the top of the lowest ice lid
        top_ice_lid_top = top_ice_lid_bot - find(flipud(ice_ind(1:top_ice_lid_bot)) == 1, 1, 'last') + 1;

        %if the top_ice_lid_top is not the top cell then there is either an
        %error or there is slush or water above and things are getting
        %ridiculous
        if any(ice_ind(1:top_ice_lid_bot) == 0)
            warning('There are 3 lakes and things are getting ridiculous.'); fprintf('%.2f', (tt*t_step)/24)
        end

        %these statements are repeated in the following elseif too as I'm too
        %tired to alter it, this makes logical sense, and it won't affect 
        %performance anyway

        %find the lowest water cell of the top lake
        temp_ind = min(strfind(fliplr(water_ind(1:bot_ice_lid_top)'), [0 1])); %temporary used as incorporated to top_lake_bot if criteria met
        
        top_lake_bot = bot_ice_lid_top - temp_ind;
  
        %find the top water cell of the top lake
        top_lake_top = find(water_ind == 1, 1, 'first');

        %calculate top lake index
        top_lake_ind = (top_lake_top:top_lake_bot)';

        %find the lowest water cell of the bottom lake
        bot_lake_bot = total_grid_num - min(strfind(fliplr(water_ind'), [0 1]));

        %find the top water cell of the bottomm lake
        bot_lake_top = bot_ice_lid_bot + find(water_ind(bot_ice_lid_bot:end) == 1, 1, 'first') - 1;

        %calculate bottom lake index
        bot_lake_ind = (bot_lake_top:bot_lake_bot)';

    %if there is ice and slush above the lowest lid
    elseif any(ice_ind(1:bot_ice_lid_top - 1) == 1) && any(slush_ind(1:bot_ice_lid_top) == 1)

        %this line will find the bottom of the top ice lid
        top_ice_lid_bot = min(ice_ind_sub_ud(find(diff(ice_ind_sub_ud) < - ice_lim, 2, 'first') + 1));

        %this finds the top of the lowest ice lid
        top_ice_lid_top = top_ice_lid_bot - find(flipud(ice_ind(1:top_ice_lid_bot)) == 1, 1, 'last') + 1;
        
        %if the top_ice_lid_top is not the top cell then there is either an
        %error or there is slush or water above and things are getting
        %ridiculous
        if any(ice_ind(1:top_ice_lid_top) == 0) 
            
            %is there water above the top ice lid?
            if any(water_ind(1:top_ice_lid_top) == 1)
                warning('There are 3 lakes and things are getting ridiculous.'); fprintf('%.2f', (tt*t_step)/24)
                
            %if it's only slush then combine lids and reset top lid. This
            %means that slush is incorporated into lid
            elseif any(slush_ind(1:top_ice_lid_top) == 1)
                
                %are there 'too many' slush cells within proposed lid
                %merger
                if bot_ice_lid_top - top_ice_lid_bot <= slush_in_lid_max
                    
                    %combine lids and reset top lid
                    bot_ice_lid_top = top_ice_lid_top;
                    top_ice_lid_top = [];
                    top_ice_lid_bot = [];
                end
                
            %if this code is accessed then there is only one lake and its
            %bounds and index can be calculated as above
            bot_lake_bot = total_grid_num - find(flipud(water_ind) == 1, 1, 'first') + 1;
            bot_lake_top = bot_ice_lid_bot + find(water_ind(bot_ice_lid_bot:end) == 1, 1, 'first') - 1;
            bot_lake_ind = (bot_lake_top:bot_lake_bot)';
            
            %set tracker to prevent recalculation of bot_lake_ind
            slush_track = 1;
            
            else
                error('An unacounted for scenario has arisen.')
            end
        end

    %if there is water or slush above the lowest lid
    elseif any(water_ind(1:bot_ice_lid_top - 1) == 1) || any(slush_ind(1:bot_ice_lid_top - 1) == 1)

        %find the lowest water cell of the top lake
        top_lake_bot = bot_ice_lid_top - min(strfind(fliplr(water_ind(1:bot_ice_lid_top)'), [0 1]));

        %find the top water cell of the top lake
        top_lake_top = find(water_ind(1:bot_ice_lid_top) == 1, 1, 'first');

        %calculate top lake index
        top_lake_ind = (top_lake_top:top_lake_bot)';

        %if top_lake_top is empty then only slush is above top ice lid. In this
        %case set just this one slush cell as lake
        if isempty(top_lake_top)
            top_lake_top = 1;
            top_lake_ind = 1;
        end

        %find the lowest water cell of the bottom lake
        bot_lake_bot = total_grid_num - min(strfind(fliplr(water_ind'), [0 1]));

        %find the top water cell of the bottomm lake
        bot_lake_top = bot_ice_lid_bot + find(water_ind(bot_ice_lid_bot:end) == 1, 1, 'first') - 1;

        %calculate bottom lake index
        bot_lake_ind = (bot_lake_top:bot_lake_bot)';

        if ~ismember(1,water_ind(bot_ice_lid_bot:ice_substrate_top))
            %warning('Only slush between lower ice lid and ice substrate'); fprintf('%.2f', (tt*t_step)/24)

        end

    end
    
    %if there is slush but no ice between the bottom of the bottom lake and
    %the surface there then must be a slush lid. Here using 2 slush cells
    %to prevent this being accessed if there is just one surface slush cell
    if sum(slush_ind(1:bot_lake_bot)) >= 2 && ~ismember(1, ice_ind(1:ice_substrate_top - 1))
        slush_lid_bot = bot_lake_bot - find(flipud(slush_ind(1:bot_lake_bot)) == 1, 1, 'first') + 1;
        slush_lid_top = find(slush_ind == 1, 1, 'first');
    
    %slush lid also needs setting as ice incase there is ice when lid
    %breakup occurs (which is the only function slush_lid_top/bot are used
    %for
    elseif ~ismember(1, ice_ind(1:ice_substrate_top - 1))
        slush_lid_bot = bot_ice_lid_bot;
        slush_lid_top = bot_ice_lid_top;
    end
    
    %IF LID IS COMPRISED ENTIRELY OF SLUSH BUT THERE ARE TWO DISTINCT LAKES
    
     %this is all a bit complicated but I think it works. Possible area
    %where errors may occur
    if sum(lambda(bot_lake_ind) <= slush_lid_threshold) >= slush_lid_num && slush_track == 0
        
        top_lake_ind = (1:find(lambda < 1, 1, 'first') - 1)';
        
        bot_lake_bot = total_grid_num - find(flipud(water_ind) == 1, 1, 'first') + 1;
        
        %find where the top slush layer ceases and bottom layer starts
        bot_lake_top = max(top_lake_ind) + find(slush_ind(max(top_lake_ind) + 1:end) == 0, 1, 'first');
        %bot_lake_top = bot_ice_lid_bot + find(water_ind(lid_assumption:end) == 1, 1, 'first') - 1
        bot_lake_ind = (bot_lake_top:bot_lake_bot)';
        
    end
    
    %calculate lake depths including bounding slush cells if present
    if ~isempty(bot_lake_ind) && ~isempty(top_lake_ind)

        %bottom lake if top lake exists
        bot_lake_depth = sum(grid_profile(bot_lake_ind));% +...
        %lambda(max(bot_lake_ind) + 1)*grid_profile(max(bot_lake_ind) + 1) +...
        %lambda(min(bot_lake_ind) - 1)*grid_profile(min(bot_lake_ind) - 1);

    elseif ~isempty(bot_lake_ind) && isempty(top_lake_ind)

        %bottom lake if top lake does not exist (and therfore no lid)
        bot_lake_depth = sum(grid_profile(bot_lake_ind)) +...
            lambda(max(bot_lake_ind) + 1)*grid_profile(max(bot_lake_ind) + 1);

    else
        bot_lake_depth = 0;
    end

    if ~isempty(top_lake_ind)

        %top lake
        top_lake_depth = sum(grid_profile(top_lake_ind));

    else
        top_lake_depth = 0;
    end


end
