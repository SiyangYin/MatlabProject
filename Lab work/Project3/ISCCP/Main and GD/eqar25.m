
%c
%c=====================================================================  
%c     DLON AND DLAT CAN BE ANY VALUES SPECIFIED IN MAIN PROGRAM.        
%c---    NOTE: DIMENSION: 1:#                                            
%      SUBROUTINE eqar25 ( DLON,DLAT)
%c                                                                       
%c     EQUAL-AREA COMPUTATIONS                                           
%c                                                                       
%**                                                                    **
%      REAL     DLAT, DLON                                               
%      INTEGER  latmax, LONMAX       

%c     PARAMETER ( DLON=0.5, DLAT=0.5 )                                  
%      PARAMETER ( LONMIN=1, latmin=1 )
      function [dlontb,ncells,loneqr,lateqr,xlatb,xlate,boxnmb]=eqar25(dlon,dlat)
      lonmin=1; latmin=1;
%      dlon=2.5; dlat=2.5;
%c     PARAMETER ( LONMAX=360.0/DLON, latmax=180.0/DLAT )                
%**                                                                    **
%      INTEGER*2 boxnmb,loneqr,lateqr,ncells     
%      common/eqreqa/ boxnmb(144,72),loneqr(6596),lateqr(6596),     
%     &               ncells(72),dlontb(72),xlatb(72),xlate(72)
        boxnmb=ones(144,72);loneqr=ones(1,6596); lateqr=ones(1,6596);
        ncells=ones(1,72); dlontb=ones(1,72); xlatb=ones(1,72); xlate=ones(1,72);
%**                                                                    **
%      DATA re /6371.2/      
      re=6371.2;
%c                                                                       
%c     BEGIN                                                             
      lonmax=360/dlon;                                                 
      latmax=180/dlat;                                                 
%c                                                                       
%     PI=2.*ASIN(1.);                                                    
      twopi=2*pi;                                                       
      twopir=twopi/360;                                                 
      rdlat=dlat*twopir;                                                 
      hezon=re*sin(rdlat);                                               
      aezon=twopi*re*hezon;                                              
%co    AECELL=(AEZON*DLAT)/360. 
      aecell=(aezon*dlon)/360;                                          
%c     PRINT 10,DLAT,latmax,AECELL                                       
%   10 FORMAT(/,2X,'DLAT',F8.2,2X,'latmax',I5,2X,'AREA EQ CELL',         
%     1  F8.1,/) 
%     fprintf(1,'  dlat %8.2f  latmax %5i  area eq cell %8.1f \n',[dlat,latmax,aecell]);
%c         
   xlatb=ones(1,latmax);   
   xlate=ones(1,latmax);
   ncells=ones(1,latmax);
      for  lat=latmin:latmax              %DO 100                           
           xlatb(lat)=dlat*(lat-1)-90.0;
           xlate(lat)=xlatb(lat)+dlat;                              
           if xlate(lat)>90
               xlate(lat)=90;
           end
           if xlatb(lat)>90
               xlatb(lat)=n:latmax;              %DO 100                           
               xlatb(lat)=dlat*(lat-1)-90.0;
           xlate(lat)=xlatb(lat)+dlat;
           end
           if xlate(lat)>90
               xlate(lat)=90;
           end
           if xlatb(lat)>90
               xlatb(lat)=90;         
           end
      rlatb=twopir*xlatb(lat);                     
      rlate=twopir*xlate(lat);                                        
      htb=re*sin(rlatb);                                        
      hte=re*sin(rlate);                                                
      htzone=abs(hte-htb);                                               
      azone=twopi*re*htzone;                                             
      cells=azone/aecell;                                             
      ncells(lat)=floor(cells+0.5);
   continue;                              %100 CONTINUE  
      end
%c                                                                       
      latlft=latmin;                                          
      if  mod(latmax-latmin+1,2)==1 
           latrit=latmax-1;                                            
      else                                                              
           latrit=latmax;
      end                                                             
      for  lat=latmin:latmax         %DO 200                                  
           if latlft>=latrit     %GOTO 201  
               continue;             %201 CONTINUE
           else
           ncells(latrit)=ncells(latlft);                            
           xlatb(latrit)=-xlate(latlft);                            
           xlate(latrit)=-xlatb(latlft);
           latlft=latlft+1;                                              
           latrit=latrit-1;     
           end
   continue;                        %200 CONTINUE
      end                                                        
%c                                                                       
%c     WRITE(31,3000)                                                    
% 3000 FORMAT(1X,'lat',5X,'xlatb',5X,'xlate',5X,'ncells',5X,'dlontb')   
%  dlontb=ones(latmax);
%      fprintf(fid_31,'lat     xlatb     xlate     ncells     dlontb',[]);
      for  lat=latmin:latmax        %DO 300                                 
      if  ncells(lat)>0                                  
           dlontb(lat)=360/ncells(lat);
           else                                                              
           dlontb(lat)=360.0;                                          
      end                                                             
      if xlate(lat)>90
          xlate(lat)=90;
%c     WRITE(21,301) lat,xlatb(lat),xlate(lat),                      
%c    &          ncells(lat),dlontb(lat)                             
%  301 FORMAT(1X,'lat',I4,1X,'XLAT',2F8.1,2X,'CELLS',                    
%     1  I8,2X,'DLON',F8.1)    
      fprintf(fid_21,'lat %4i xlat %8.1f %8.1f  cells %8i  dlon %8.1f \n',[lat,xlatb(lat),xlate(lat),ncells(floor(lat)),dlontb(floor(lat))]);
      end
   continue;                    %300 CONTINUE
      end
    
%c::   Added from cveq25()
   boxnmb=ones(144,72);
   loneqr=ones;
   lateqr=ones;
      for  j=1:72              %DO 401                                        
      for  i=1:144             %DO 401                                       
          boxnmb(i,j)=-1000;                                             
   continue;                   %401CONTINUE           
      end
      end
          KOUNT=0;                                                       
      for  j=1:72                %DO 411                                     
      for  i=1:ncells(j)         %DO 411                                     
          KOUNT=KOUNT+1;                                                 
          boxnmb(i,j)=KOUNT;                                             
          loneqr(KOUNT)=i;                                               
          lateqr(KOUNT)=j;                                               
   continue;                     %411 CONTINUE     
      end
      end
%CC    IF(KOUNT .NE. 6596)STOP 4444  
%c::
      return;      