function [MAP std_dev] = plot_posterior( posterior_sample, N_bins, variable_name, period )
%PLOT_POSTERIOR makes a histogram of posterior
% returns maximum a posteriori

% posterior_sample      posterior sample
% N_bins                number of in histogram
% variable_name         variable name in plot title
% period                period of variable (use 0 for non-periodic)

% MAP                   maximum a posteriori
% std_dev               1 sigma one-sided error of MAP

[bin_count bin_center] = hist(posterior_sample, N_bins);

[bin_count_sorted indices] = sort(bin_count);
MAP = bin_center(indices(end));

sorted_diff_from_MAP = sort( abs(posterior_sample - MAP) );
sorted_diff_from_MAP = sorted_diff_from_MAP - ( sorted_diff_from_MAP > period ) * period;
sorted_diff_from_MAP = sort(sorted_diff_from_MAP);
std_dev = sorted_diff_from_MAP( ceil(length(posterior_sample)*0.68) );

bar(bin_center,bin_count,1, 'FaceColor', 'none', 'EdgeColor', 'b','linewidth',1)
hold on
line([MAP MAP],[0 2*max(bin_count)],'color','r','linewidth',2)
hold off
bin_width = bin_center(2)-bin_center(1);
axis([min(bin_center)-bin_width max(bin_center)+bin_width 0 1.2*max(bin_count)]);
title_str = ['$' variable_name '=' num2str(MAP,4) '\pm' num2str(std_dev,2) '$'];
title(title_str,'interpreter','latex','fontsize',10)
set(gca,'ytick',[], 'yticklabel',{})

end

