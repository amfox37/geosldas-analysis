function [ stats, stats_tags ] = get_validation_stats( data, AC, complete_or_pairwise, ...
                        ref_col, select_col, N_pairs_min, get_CI ) 

% calculate validation statistics for the incoming data-matrix 
%
% input:  data :       [NTime, (insitu SMOS exp1 exp2...)] 
%         AC:          0, 1 (autocorrelated=1 or not)
%         complete_or_pairwise 
%         ref_col:     reference (insitu) data for validation
%         select_col:  only consider data prior to masking columns
%         N_pairs_min: minimum # of data pairs
%         get_CI:      0 or 1 (0 would avoid looking for unavailable
%         statistics toolbox)
%
% output: stats
%         stats_tags
%
% GDL, 22 Jan 2014: re-organized
% GDL, 09 Oct 2015: added CI for ubMSE
% GDL, 13 Oct 2015: bug-fixes, 
%                   consistently deal with adjusted sample sizes (autocorrelation)
%                   inside the subroutines with statistical tests 
%
% Documentation on CI:
% http://www.dtcenter.org/met/users/docs/users_guide/MET_Users_Guide_v4.1.pdf
% ---------------------------------------------------------------------

ref_col = find(ref_col == select_col);

if isempty(ref_col) || ref_col~=1
  error('reference (insitu) data need to be the first entry of the selected data columns')
end

nodata_val = -9999;
tol        = 1E-04;

%stats_tags = {'N_pairs','R', 'RLO', 'RUP','bias','CI_bias','MSE', 'CI_MSE'};
%changed 9 dec 2014
%stats_tags = {'N_pairs','R', 'RLO', 'RUP','bias','CI_bias','MSE', 'CI_MSE','ubMSE'};
%changed 9 oct 2015
stats_tags = {'N_pairs','R', 'RLO', 'RUP','bias','CI_bias','MSE', 'CI_MSE','ubMSE', 'ubMSELO', 'ubMSEUP'};
N_stats    = length(stats_tags);
N_select   = length(select_col)-1; %all stats are versus the ref_col

stats      = nodata_val+zeros(N_select,N_stats);
stats(:,1) = 0 ;

%----------------------------------------------------------------------

%% 1 = Number of pairs 
%   Reduce data to exact overlapping sets.
%   Here the masking versus the later excluded columns is included.

if (strcmp(complete_or_pairwise,'complete'))

    t = isnan(data);
    data = data(~any(t,2),:);

else %typically when the ass-obs are not avail

    t1 = isnan(data(:,ref_col));       %insitu
    t2 = isnan(data(:,select_col(3))); %model - dangerously hardwired; select_col could be anything!
    t = t1;
    t((t2==1)) = 1;
    data = data(~any(t,2),:);

end

stats(:,1) = size(data,1);

%cut back columns that were only meant for crossmasking
data = data(:,select_col);

%----------------------------------------------------------------------

if size(data,1)>N_pairs_min

  %% 2,3,4 = Correlation

  %If a full time series is used for validation: pay attention to
  %autocorrelation!
        %To get the effective sample size for the confidence limit
        %calculation, we need to account for the autocorrelation in the
        %datasets; since the datasets contain NAN, we need a clever way to
        %get the lag-1 correlation, i.e. through rescaling and setting NaNs
        %to zero. => all incorporated into corrcoef_autocorr

  %The confidence bounds are based on an asymptotic normal distribution of 0.5*log((1+R)/(1-R)), with an approximate variance equal to 1/(n-3). These bounds are accurate for large samples when X has a multivariate normal distribution.

   if (AC)
         [Rtmp,Ptmp,RLOtmp,RUPtmp]= corrcoef_autocorr(data, 'rows', complete_or_pairwise);
   else
         [Rtmp,Ptmp,RLOtmp,RUPtmp]= corrcoef(data, 'rows', complete_or_pairwise);
   end

   if (size(Rtmp) > 1)

        %% 2 = Correlation
        stats(:,2)  = Rtmp(1,select_col(2:end));

        %% 3 = Lower bound R
        stats(:,3)  = RLOtmp(1,select_col(2:end)); 

        %% 4 = Upper bound R
        stats(:,4)  = RUPtmp(1,select_col(2:end));

    else

        stats(:,2:4)= nodata_val;

    end


   %------------------------------------------------------------------
   %  Calculate sample_size for CI of bias and RMSE, based on individual
   %  x- and y-time series (Neff = [1-rxry)(1+rxry)])
   %  to avoid a too small effective sample when
   %  calculating Neff based on the differences (which are highly
   %  autocorrelated)
   %  (Neff = [1-rdiff^2)(1+rdiff^2)]
   %  - ceil the sample size to moderate the very small sample sizes
   %    for highly correlated time series (as for 3-hourly data)
   %
   if (AC)
       sample_size = zeros(1,length(select_col)-1);

       for i=1:length(select_col)-1
        xx = data(:,ref_col);
        yy = data(:,select_col(ref_col+i));

        sample_size(1,i) = sum(~isnan(xx) & ~isnan(yy));

        rx = nancorrcoef(xx(1:end-1),xx(2:end));
        ry = nancorrcoef(yy(1:end-1),yy(2:end));

        sample_size(1,i) = ceil(sample_size(1,i).*(1-rx*rx)./(1+rx*rx));
       end
   else
       sample_size = NaN; %b/c detected inside statistical test
   end
   
   %------------------------------------------------------------------
   %% 5,6 = Bias

   err = (repmat(data(:,ref_col),1,N_select) -...
                   data(:,select_col(2:end)) ); %(insitu-model)
   stats(:,5) = nanmean( err, 1);

   if get_CI
       %function adjusted for:
       % - sample size in T-test
       % - deal with non-zero mean which is not a scalar across all
       %   experiments: this is only important if p-values or other stuff 
       %   is requested; otherwise the '0' will return just the same CI
       [h_T,p_T,ci] = Ttest_signf_autocorr(err,stats(:,5)','Alpha',0.05,'sample_s',sample_size);
   else
       ci = NaN*ones(2,size(stats,1));
   end

   %ci(:,1) = lower bound, ci(:,2) = upper bound (equal distance from mean)
   %ci(2,:) = upper bound
   stats(:,6) = ci(2,:)-nanmean(err, 1);

   %------------------------------------------------------------------
   %% 7,8 = MSE

   err2 = err.^2;
   stats(:,7) = nanmean( err2, 1);
 
   if get_CI
       %null hypothesis: mean = stats(:,7)
       [h_T,p_T,ci] = Ttest_signf_autocorr(err2, stats(:,7)' ,'Alpha',0.05,'sample_s',sample_size);
   else
       ci = NaN*ones(2,size(stats,1));
   end
  
   %CI have equal distance added/substracted to/from mean
   %ci(1,:) = lower bound 
   %ci(2,:) = upper bound
   stats(:,8) = ci(2,:)-stats(:,7)';   
 
   %--------------------------------------------------------------------
   %% 9 = ubMSE = MSE - bias^2
   %% 9,10,11 = MS unbiased E (=variance)

   stats(:,9) = stats(:,7) - stats(:,5).^2;
   err_0      = err - repmat(stats(:,5)',size(err,1),1);
   
   if abs(squeeze(stats(:,9)') - nanmean( err_0.^2, 1)) > tol
       disp('something is wrong with unbiased MSE');
   end
  
   for i=1:size(err2,2)
       if get_CI
           [h_T,p_T,ci] = vartest_autocorr(err_0(:,i), stats(i,9)','Alpha',0.05,'sample_s',sample_size);
        else
           ci = NaN*ones(2,size(stats(:,i),1));
       end
       
       %CI have NON-equal distances around variance
       %ci(1,:) = lower bound 
       %ci(2,:) = upper bound 
       stats(i,10) = ci(1);
       stats(i,11) = ci(2);

  end
    
   %------------------------------------------------------------------
   %limit to T-test output for MSE, add ubMSE
   %stats = stats(:,1:8);
   %stats = stats(:,1:9);

end


% ========= EOF =========================================================
  
