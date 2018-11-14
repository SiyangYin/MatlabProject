%{
c=====================================================================  
*-------------------------------------------------------------------    
*           REPLICATING REAL*4 ARRAY (6596) TO (144,72) STARTING        
*     AT EITHER GREENWICH LINE (ioptn=0) OR DATELINE (ioptn=1)!         
*-------------------------------------------------------------------    
%}
%{      
                SUBROUTINE eq25ea(ioptn,arybox,eqangl)                 
      REAL arybox(6596),eqangl(144,72)                                  
      INTEGER*2 boxnmb,loneqr,lateqr,ncells     
      common/eqreqa/ boxnmb(144,72),loneqr(6596),lateqr(6596),     
     &               ncells(72),dlontb(72),xlatb(72),xlate(72)
 %}
%c--- 
      function [boxnmb,eqangl]= eq25ea(ioptn,arybox,eqangl)
%      ioptn=0; arybox=zeros(1,7000); arybox(1:25)=[2409,2151.20000000000,0.639600000000000,-1000,13.5680000000000,0,-1000,3098.50000000000,11.0080000000000,23.8100000000000,972.800000000000,-1000,-1000,0,0,1664.50000000000,-1000,0,0,-1000,-1000,0,0,-1000,-1000];
      [dlontb,ncells,loneqr,lateqr,xlatb,xlate,boxnmb]=eqar25(2.5,2.5);    
      for  j=1:72           %DO 101                                          
      for  i=1:144          %DO 101                                          
          eqangl(i,j)=-1000;                                            
                            %101 CONTINUE    
      end
      end
      for  j=1:72             %DO 111                                        
          lat=j;                                                        
      for  i=1:144            %DO 111                                        
          if(ioptn==0)                                           
             iadj=i;                                                     
          else if ioptn==1                                     
             iadj=i+72; 
              end
          end
             if(iadj>144)
                 iadj=iadj-144;
             end                                                        
          lon=floor(((i-0.5)*2.5)/dlontb(lat)) + 1;
          if lon>ncells(lat)
          lon=ncells(lat);     
          end
          eqangl(iadj,j)=arybox(boxnmb(lon,j));
      end
      end
                            %111 CONTINUE
      
       return;     
      end