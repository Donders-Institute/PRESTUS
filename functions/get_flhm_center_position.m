function [flhm_center, flhm_center_index] = get_flhm_center_position(x,y)
  halfMax = (min(y) + max(y)) / 2;
  [~,max_pos] = max(y);
  % Find where the data first drops below half the max.
  index1 = find(y <= halfMax & x < x(max_pos), 1, 'last')+1;
  % Find where the data last rises above half the max.
  index2 = find(y <= halfMax & x > x(max_pos), 1, 'first')-1;

  flhm_center = (x(index2) - x(index1))/2+x(index1);

  [~, flhm_center_index] = min(abs(x-flhm_center));
end