function date = Date2Julian(year, days, opt)
    arguments
        year;
        days;
        opt.month(1,1) = 0;
        opt.hour (1,1) = 0;
        opt.min (1,1) = 0;
        opt.sec (1,1) = 0;
    end

    if  opt.month == 0
        switch days
            case days*((days>0) & (32 > days))
               opt.month= 1;
               day = floor(days);
            case days*((days>=32) & (60 > days))
               opt.month= 2;
               day = floor(days-31);
            case days*((days>=60) & (91 > days))
               opt.month= 3;
               day = floor(days-59);
            case days*((days>=91) & (121 > days))
               opt.month= 4;
               day = floor(days-90);
            case days*((days>=121) & (152 > days))
               opt.month= 5;
               day = floor(days-120);
            case days*((days>=152) & (182 > days))
               opt.month= 6;
               day = floor(days-151);
            case days*((days>=182) & (213 > days))
               opt.month= 7;
               day = floor(days-182);
            case days*((days>=213) & (244 > days))
               opt.month= 8;
               day = floor(days-212);
            case days*((days>=244) & (274 > days))
               opt.month= 9;
               day = floor(days-243);
            case days*((days>=274) & (304 > days))
               opt.month= 10;
               day = floor(days-273);
            case days*((days>=304) & (335 > days))
               opt.month= 11;
               day = floor(days-303);
            case days*((days>=335) & (366 > days))
               opt.month= 12;
               day = floor(days-334);
        end
        time_rem = rem(days,1)*86400;
        opt.hour = floor(time_rem/3600);
        opt.min = floor(time_rem/60)-opt.hour*60;
        opt.sec = time_rem-opt.hour*60^2-opt.min*60;
        days = day;
    end

    comp1 = 367*year;
    comp2 = floor(7*(year+floor((opt.month+9)/12))/4);
    comp3 = floor(275*opt.month/9);
    comp4 = days;
    comp5 = (((opt.sec/60+opt.min)/60)+opt.hour)/24;
    date = comp1 - comp2 + comp3 + comp4 + 1721013.5 + comp5;
end