function STA = sim_sta(x)
%{
    Compute differences to the first sample and basic statistics.
    INPUT:
        x : Repeating arguments (arrays/matrices), interpreted as {x1, x2, x3, ...}
    OUTPUT:
        STA : struct with fields
            .dif    : column-stacked differences [vec(x2-x1), vec(x3-x1), ...]
            .no_nan : number of non-NaN entries per column (same size as columns in dif)
            .max    : per-column maximum (NaNs propagate)
            .min    : per-column minimum (NaNs propagate)
            .mean   : per-column mean (ignoring NaNs)
            .std    : per-column std (ignoring NaNs)
            .meadad : per-column median absolute deviation (MAD, about median, flag=1)
            .meanad : per-column mean absolute deviation (MAD default)
            .mae    : mean absolute error (all non-NaN dif values concatenated)
            .mare   : mean absolute relative error (%) w.r.t x1 (NaN-masked)
            .rms    : per-column RMS (normalized by non-NaN count in that column)
            .rrms   : relative RMS (column RMS divided by RMS(x1))
            .coff   : Pearson correlation between x1 and (dif + x1) per column
%}
arguments(Repeating)
    x
end

num = length(x);
x1 = x{1};
len = length(x1(:));

dif = zeros(len, num - 1);

if num > 1
    for i = 1 : num - 1
        dif(:, i) = reshape(x{i + 1} - x1, [], 1);
    end
else
    dif = x1(:);
end

STA = struct();
STA.dif = dif;

ind = ~isnan(dif);
STA.no_nan = sum(ind);
STA.mean = nanmean(dif);
STA.std = nanstd(dif);
STA.max = nanmax(dif);
STA.min = nanmin(dif);
STA.meadad= mad(dif, 1);
STA.meanad= mad(dif);
dif(~ind) = [];
STA.mae = mean(abs(dif));
x1 = x1(:);
x1(~ind) = [];
STA.mare = mean(abs(dif./x1)) * 100;
STA.rms = sqrt(sum(dif .^ 2) ./ STA.no_nan);
STA.rrms = STA.rms / rms(x1);

if num == 1
    STA.coff = [];
else
    STA.coff = corr(x1, dif + x1, 'type', 'Pearson');
end

end

