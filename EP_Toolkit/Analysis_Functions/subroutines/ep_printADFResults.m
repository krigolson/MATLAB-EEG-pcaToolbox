function ep_printADFResults(RESULTS, alpha, followup, outfid, NUMSIM);
%ep_printADFResults(RESULTS, alpha, followup, outfid, NUMSIM);
%Prints out statistical results with appropriate formatting.
%
%Inputs
%  alpha :   The alpha thresholds for significance
%    .uncorrected:  Alpha for a priori interest (.05)
%    .corrected:  Alpha for posthoc multiple comparison corrected.
%    .pvar:  twice standard deviation of p-values over reps
%  RESULTS	: (1) = test statistic
%  RESULTS	: (2) = Numerator DF
%  RESULTS	: (3) = Denominator DF
%  RESULTS	: (4) = Significance
%  RESULTS	: (5) = Effect size using delta-hat-star standardizer.
%  RESULTS	: (6) = Lower confidence limit
%  RESULTS	: (7) = Upper confidence limit
%  RESULTS	: (8) = Scaling factor for effect size estimator if SCALE option is chosen.
%  RESULTS	: (9) = Number of singular matrices during calculation of statistics.
%  RESULTS	: (10) = Number of nearly singular matrices during calculation of statistics.
%  outfid:  The fid of the output file.  Will print to screen if zero.
%  NUMSIM   : Number of simulation runs in the bootstrapping procedure.
%Outputs
%  Prints out statistical results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%History
%  by Joseph Dien (11/8/08)
%  jdien07@mac.com
%
% modified 2/12/09 JD
% Correctly handles NaN results.
%
% modified 3/23/09 JD
% Highlights significant results.  Added option to print to file.
%
% modified 4/21/09 JD
% Green now corresponds to significant 1-tailed test (two times the uncorrected alpha).
%
% modified 7/26/09 JD
% Added indentation for follow-up ANOVAs.
%
% modified 9/20/15 JD
% Added effect sizes.
%
% bugfix 11/30/15 JD
% No longer produces effect sizes for contrasts of more than 1 df.
%
% modified 12/11/15 JD
% Checks if p-value variability exceeds two standard deviations over the threshold.
% p-values exactly equalling the threshold are accepted as significant.
% Prints out number of singular and nearly singular matrices during computation of statistics.
%
% modified 12/27/15 JD
% Changed message to indicate that effect size only available for
% between-group contrasts.

if nargin < 4
    outfid =0;
end;

P=sprintf('%10.8f', RESULTS(4));
Pp=findstr('.',P);
Pz=findstr(P(Pp+1:length(P)),'0');
for i=1:length(Pz)
    if i ~= Pz(i)
        break
    end
end
if isempty(i)
    i=1;
end;

if RESULTS(4) >= .00000001
    pSign='=';
    pVal=sprintf(['%1.' num2str(i+1) 'f'], RESULTS(4));
elseif isnan(RESULTS(4))
    pSign='=';
    pVal=' could not be calculated';
else
    pSign='<';
    pVal='0.00000001';
end

if RESULTS(4) <= alpha.corrected
    if (RESULTS(4)+alpha.pvar) > alpha.corrected
        pVal=[pVal '(not confirmed)'];
    end;
elseif RESULTS(4) <= alpha.uncorrected
    if (RESULTS(4)+alpha.pvar) > alpha.uncorrected
        pVal=[pVal '(not confirmed)'];
    end;
end;

if isnan(RESULTS(5))
    effsz=sprintf(', effect sizes only for between-group one degree of freedom contrasts'); 
else
    effsz=sprintf(', d*=%1.1f, CI<sub>95</sub>=%1.2f to %1.2f', RESULTS(5), RESULTS(6), RESULTS(7)); 
end;


if outfid
    %print to file
    if RESULTS(4) <= alpha.corrected
        fprintf(outfid,[repmat('&nbsp;',1,followup*2) '<font color="red">T<sub>WJt</sub>/c(%1.1f,%1.1f)=%1.2f, p%s%s%s.</font></BR>\n'], RESULTS(2), RESULTS(3), RESULTS(1), pSign, pVal, effsz);
    elseif RESULTS(4) <= alpha.uncorrected
        fprintf(outfid,[repmat('&nbsp;',1,followup*2) '<font color="orange">T<sub>WJt</sub>/c(%1.1f,%1.1f)=%1.2f, p%s%s%s.</font></BR>\n'], RESULTS(2), RESULTS(3), RESULTS(1), pSign, pVal, effsz);
    elseif RESULTS(4) <= 2*alpha.uncorrected
        fprintf(outfid,[repmat('&nbsp;',1,followup*2) '<font color="green">T<sub>WJt</sub>/c(%1.1f,%1.1f)=%1.2f, p%s%s%s.</font></BR>\n'], RESULTS(2), RESULTS(3), RESULTS(1), pSign, pVal, effsz);
    else
        fprintf(outfid,[repmat('&nbsp;',1,followup*2) 'T<sub>WJt</sub>/c(%1.1f,%1.1f)=%1.2f, p%s%s%s.</BR>\n'], RESULTS(2), RESULTS(3), RESULTS(1), pSign, pVal, effsz);
    end;
    if RESULTS(9) > 0
        fprintf(outfid,'Warning: The matrix was singular to working precision %4.2f%%%s',100*(RESULTS(9)/NUMSIM), ' of the time when computing the significance.  These runs were ignored and thus the estimate may have some inaccuracy.</BR>');
    end;
    if RESULTS(10) > 0
        fprintf(outfid,'Warning: Matrix was singular, close to singular or badly scaled %4.2f%%%s',100*(RESULTS(10)/NUMSIM), ' of the time when computing the significance.  These runs were ignored and thus the estimate may have some inaccuracy.</BR>');
    end;
else
    %print to screen
    if RESULTS(4) <= alpha.corrected
        fprintf([repmat('&nbsp;',1,followup*2) '**TWJt/c(%1.1f,%1.1f)=%1.2f, p%s%s%s.**\n'], RESULTS(2), RESULTS(3), RESULTS(1), pSign, pVal, effsz);
    elseif RESULTS(4) <= alpha.uncorrected
        fprintf([repmat('&nbsp;',1,followup*2) '*TWJt/c(%1.1f,%1.1f)=%1.2f, p%s%s%s.*\n'], RESULTS(2), RESULTS(3), RESULTS(1), pSign, pVal, effsz);
    else
        fprintf([repmat('&nbsp;',1,followup*2) 'TWJt/c(%1.1f,%1.1f)=%1.2f, p%s%s%s.\n'], RESULTS(2), RESULTS(3), RESULTS(1), pSign, pVal, effsz);
    end; 
    if RESULTS(9) > 0
        fprintf(['Warning: The matrix was singular to working precision ' sprintf('%4.2f%%',100*(RESULTS(9)/NUMSIM)) ' of the time when computing the significance.  These runs were ignored and may thus the estimate may have some inaccuracy.\n']);
    end;
    if RESULTS(10) > 0
        fprintf(['Warning: Matrix was singular, close to singular or badly scaled ' sprintf('%4.2f%%',100*(RESULTS(10)/NUMSIM)) ' of the time when computing the significance.  These runs were ignored and may thus the estimate may have some inaccuracy.\n']);
    end;
end;
