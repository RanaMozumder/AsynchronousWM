function [AlignedDelayRates, AlignedClasses, AlignedDelaySEM, SpkCount] = SpikeCountVsClass(filename, epoch, varargin)
%%
load(filename)
switch numel(varargin)
    case 0
        trials = 'all';
    case 1
        trials = varargin{1};
    otherwise
        error('Too many arguments.')
end
if ~isempty(MatData)
    % calculate best class
    if length(MatData.class) == 8 || length(MatData.class) == 120
        for c = 1:length(MatData.class)
            if ~isempty(MatData.class(c).ntr)
                ntr = 0;
                for tr = 1:length(MatData.class(c).ntr)
                    if ismember(MatData.class(c).ntr(tr).trialnum, trials) || strcmpi(trials, 'all')

                        ntr = ntr + 1;

                        TS       = MatData.class(c).ntr(tr).TS;
                        Fix_inT  = MatData.class(c).ntr(tr).Fix_inT;
                        Cue_onT  = MatData.class(c).ntr(tr).Cue_onT;
                        Cue_offT = MatData.class(c).ntr(tr).Cue_offT;
                        Fix_offT = MatData.class(c).ntr(tr).Fix_offT;

                        fixrate       = length(find((TS>Fix_inT)  & (TS<=Cue_onT)))/(Cue_onT-Fix_inT);
                        fixcount      = length(find((TS>Fix_inT)  & (TS<=Cue_onT)));

                        cuerate       = length(find((TS>Cue_onT)  & (TS<=Cue_offT)))/(Cue_offT-Cue_onT);
                        cuecount      = length(find((TS>Cue_onT)  & (TS<=Cue_offT)));

                        cuedelayrate  = length(find((TS>Cue_offT) & (TS<=Fix_offT)))/(Fix_offT-Cue_offT);
                        cuedelaycount = length(find((TS>Cue_offT) & (TS<=Fix_offT)));

                        rates(ntr, 1) = fixrate;
                        rates(ntr, 2) = cuerate;
                        rates(ntr, 3) = cuedelayrate;
    
                        counts(ntr, 1) = fixcount;
                        counts(ntr, 2) = cuecount;
                        counts(ntr, 3) = cuedelaycount;
    
                        clear TS fixrate cuerate cuedelayrate fixcount cuecount cuedelaycount;
                    end
                end
                Allrates(c, 1) = {rates};
                Allrates(c, 2) = {mean(rates(:,1), 'omitnan')}; % fixrate
                Allrates(c, 3) = {mean(rates(:,2), 'omitnan')}; % cuerate
                Allrates(c, 4) = {mean(rates(:,3), 'omitnan')}; % cuedelayrate
                Allrates(c, 5) = {std(rates(:,3), 'omitnan')./sqrt(mean(sum(~isnan(rates(:,3)))))}; % cuedelayrate sem

                Allcounts(c, 1) = {counts};
                Allcounts(c, 2) = {mean(counts(:,1), 'omitnan')}; % fixcount
                Allcounts(c, 3) = {mean(counts(:,2), 'omitnan')}; % cuecount
                Allcounts(c, 4) = {mean(counts(:,3), 'omitnan')}; % cuedelaycount

                clear rates counts
            
            else
                Allrates(c, 1) = {nan};
                Allrates(c, 2) = {nan}; % fixrate
                Allrates(c, 3) = {nan}; % cuerate
                Allrates(c, 4) = {nan}; % cuedelayrate
                Allrates(c, 5) = {nan}; % cuedelayrate sem
    
                Allcounts(c, 1) = {nan};
                Allcounts(c, 2) = {nan}; % fixcount
                Allcounts(c, 3) = {nan}; % cuecount
                Allcounts(c, 4) = {nan}; % cuedelaycount
            end
        end
    else
        disp('All classes unavailble')
    end
%     % best cueclass
%     [max_cuerate, max_cueclass] = max([Allrates{:,3}]); 
%     circCuemat = [Allrates{:,3}, Allrates{:,3}, Allrates{:,3}];
%     AlignedCueRates = circCuemat(max_cueclass+8-4: max_cueclass+8+4);


    % best cuedelay class
%     [~, max_cuedelayclass] = max([Allrates{:,4}]); 

%     circDelaymat = [Allrates{:,4}, Allrates{:,4}, Allrates{:,4}];
%     AlignedDelayRates = circDelaymat(max_cuedelayclass+8-4: max_cuedelayclass+8+4);
% 
%     circDelayclasses = [1:8, 1:8, 1:8];
%     AlignedClasses = circDelayclasses(max_cuedelayclass+8-4: max_cuedelayclass+8+4);
% 
%     circDelaySEM = [Allrates{:,5}, Allrates{:,5}, Allrates{:,5}];
%     AlignedDelaySEM = circDelaySEM(max_cuedelayclass+8-4: max_cuedelayclass+8+4);
    if strcmpi(epoch, 'Cue')
        SpkCount   = [Allcounts{:,3}];
    elseif strcmpi(epoch, 'Delay')
        SpkCount = [Allcounts{:,4}];
    end
% else
    AlignedDelayRates = [];
    AlignedClasses = [];
    AlignedDelaySEM = [];
%     SpkCount = [];
end