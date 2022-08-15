function settings = get_two_whisker_settings
    [irr hostname] = system('hostname');
    hostname = strtrim(hostname);
    if(length(strfind(hostname,'cluster')) == 1)
        hostname = 'cluster';
    end

    %% This dictates the location of files - hostname specific so ppl butt heads
    switch hostname
        case 'apollo'
            dat_root = '/Users/speron/Desktop/voelcker_et_al_two_whisker';
            settings.data_path_root = [dat_root '/data'];
            settings.whisking_data_path_root = [ dat_root '/whiskvid']; 
            settings.summary_data_root = [dat_root '/two_whisker_data'];

        otherwise
            error(['Your computer, ' hostname ', has no paths set; go into get_two_whisker_settings.m and make a new entry for your computer in the hostname switch statement at the start of that file']);
    end

    %% This is animal specific stuff - do not touch!
    animals = {'an279608','an280201','an283544','an284891','an284893','an274688','an279029'};
    symbols = {'^','o','v','d','s','p','h'};
    dz =       [20          20          20       20         20         20           20]; % dz between planes NOT SURE THIS IS RIGHT
    whiskers = [2 3 ;    2 3 ;      2 3 ;      2 3 ;       1 2;       1 2;        2 3]; % C row duh

    % this offset is applied to all depths because depths in data are not L1 (start of brain) normalized
    z_offset_um = [0        0          -20        -5           -40      0           0];

    % this is layer border *AFTER OFFSET* - so real 
    l2l3_um  = [171         202         186     198         201         215         204];  % MIDPOINT  between 0.1 %ile cell depth and l3l4 border
    l4l5_um  = [420         440        Inf       Inf        Inf        450         450];
    l3l4_um  = [270        310        300        310          300       320      320] ;

    subvol_offset = [0 0 0 0 0 11 0]; % most are 01 02 ... but 74688 is 12 13 ...

    for a=1:length(animals)
        settings.animals(a).name = animals{a};
        settings.animals(a).dz = dz(a);
        settings.animals(a).whiskers = whiskers(a,:);
        settings.animals(a).l2l3_border = l2l3_um(a);
        settings.animals(a).l3l4_border = l3l4_um(a);
        settings.animals(a).l4l5_border = l4l5_um(a);
        settings.animals(a).z_offset_um = z_offset_um(a);
        settings.animals(a).subvol_offset = subvol_offset(a);
        settings.animals(a).symbol = symbols{a};
    end

    settings.colors.w1Color = [1 0.1 0.1];
    settings.colors.w2Color = [0.1 0.5 1];
    settings.colors.mwColor = [1 0 1];
    settings.colors.w1w2Color = [1 0 .5];
    settings.colors.w2w1Color = [.5 0 1];
    settings.colors.proColor = [.25 .75 0.75]; % protraction color
    settings.colors.retColor = [1 0.75 .5]; % retraction color
    settings.colors.bidiColor = [0.1 0.4 .6]; % bidirectional cell color
    settings.colors.swColor = [0 0 0]; % single whisker
    settings.colors.uniColor = [0.5 0.5 .5]; % unirectional cell color

    settings.colors.mwColor = [230 128 255]/255;
    settings.colors.uniColor = [124 179 195]/255;
    settings.colors.bidiColor = [26 102 153]/255;
