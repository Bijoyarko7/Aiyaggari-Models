function idx = binarySearch(g, num)
   l = 1;
   r = length(g);
   idx = 1;
   while l < r
      idx = 1 + floor((l + r - 1) / 2);
      if g(idx) > num
        r = idx - 1; 
      elseif g(idx) <= num
        l = idx;
      end
   end
   if l == r
      idx = r; 
   end
   if g(idx) > num
     idx = -1;
   end
end