function  [h]=transparanterrorbars(x,y,yerror,color,yerror2)
% plots transparant errorbars
% yerror is half width of error bar

 if nargin ==4
     TypeError='onesided';
 elseif nargin == 5
     TypeError='doublesided';
 else
     narginchk(4,5)
 end
 % check on size
 [a,b]=size(x);
 [a1,b1]=size(y);
 [a2,b2]=size(yerror);
 
 % check on whether dimensions are similar
 if a1 == a && b1 == b
     % correct, nothing to change
 elseif a1 == b && b1 == a
     y=y';
 else
     error('x and y vectors do not have the same number of elements')
 end
 
 if a2 == a && b2 == b
     % correct, nothing to change
 elseif a2 == b && b2 == a
     yerror=yerror';
     if strcmp(TypeError,'doublesided')
        yerror2=yerrror2'; 
     end
 else
     error('x and y error vectors do not have the same number of elements')
 end
 
 %size(x)
 if (a==1) 
     
 elseif (b==1)
     x=x';
     y=y';
     yerror=yerror';
     if strcmp(TypeError,'doublesided')
         yerror2=yerror2';
     end
 else
     error('wrong size t, supply column or row vector')
 end
 
 
 if strcmp(TypeError,'onesided')
 % make polygon
        xax=[x , fliplr(x)];
        yax=[(y-yerror) , (fliplr(y+yerror))]; % first y - error , back with y y + error
 elseif strcmp(TypeError,'doublesided')
      % make polygon

        xax=[x , fliplr(x)];
        yax=[(yerror2) , (fliplr(yerror))]; % first y = error , back with y = error2 
        
 end

        % make fill
        h=fill(xax,yax,color,'FaceAlpha',0.2,'EdgeColor','none');


end

