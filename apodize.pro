PRO APODIZE,im,perc
;
;  apodization
;
     s=size(im)
     edge=100./perc
     av=mean(im)
     im=im-av
     xmask=fltarr(s(1))+1.
     ymask=fltarr(s(2))+1.
     smooth_x=rfix(s(1)/edge) ; width of the edge in x-dimension
     smooth_y=rfix(s(2)/edge) ; width of the edge in y-dimension
     print,'Now smoothing edges.'
;
;  smoothing with a cosine
;
    for i=0,smooth_x do xmask(i)=(1.-cos(!pi*float(i)/float(smooth_x)))/2.
    for i=0,smooth_y do ymask(i)=(1.-cos(!pi*float(i)/float(smooth_x)))/2.
    xmask(s(1)-smooth_x:s(1)-1)=reverse(xmask(1:smooth_x))
    ymask(s(2)-smooth_y:s(2)-1)=reverse(ymask(1:smooth_y))
    for i=0,s(1)-1 do im(i,*,*)=im(i,*,*)*xmask(i)
    for i=0,s(2)-1 do im(*,i,*)=im(*,i,*)*ymask(i)    
    im=im+av
;
    end
