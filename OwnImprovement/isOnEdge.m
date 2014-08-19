function sign = isOnEdge(mask,i,j)

if(i==1||i==size(mask,1)||j==1||j==size(mask,2))
    sign = false;
else
   region = mask(i-1:i+1,j-1:j+1);
   if(sum(region(:))>0)
       sign = true;
   else
       sign = false;
   end
end