function new_y = bounded_interp1(x,y,new_x)
 new_y = interp1(x,y,new_x,'linear',NaN);
 new_y(new_x<min(x)) = y(x==min(x));
 new_y(new_x>max(x)) = y(x==max(x));
end