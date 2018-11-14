%{
************************************ y.-c. Zhang ************ 08/2009 ****
*   Program: read_i2_FD_srf.f -- Revised from read_i2_srf.f
* to read Iintegr*2 ISCCP FD-SRF [RadFLux SRF fluxes: product (b) - see below]
**************************************************************************
*  For overall ISCCP-FD information, please read the document: i2_FD.readme 
**************************************************************************
*        The ISCCP-FD RadFlux dataset contains radiative fluxes and a 
* summary of the physical quantities used to calculate them. Shortwave and
* Longwave Radiative Flux Profiles are calculated using the following datasets 
* to specify the properties of Earth's Atmosphere and Surface: ISCCP cloud
* dataset, TOVS operational sounder products, TOMS ozone products, a 
* climatology of near-surface air temperature diurnal cycle constructed from  
* the WWW Surface Weather reports and the first NCEP reanalysis, a climatology
* of cloud particle sizes from Han et al., a climatology of stratospheric
* and upper tropospheric water vapor and stratospheric aerosols from SAGE-II, 
* a climatology of tropospheric aerosols used in the GISS climate model,
* and land surface albedo spectral dependence and spectral emissivities 
* from the GISS climate model. All of this information is collected
* into four Real*4 data products, which are reformatted to four secondary
* Integer*2 data products [(a) to (d), below], of which, a fifth monthly-mean 
* product [below (e)] is now derived from FD-PRF [(c)] (see also 
* i2_FD.readme). The following table describes the five products:
* 
* ------------------------------------------------------------------------- 
* Product Name |  Definition of Secondary Integer*2 ISCCP-FD RadFlux Products
* --------------------------------------------------------------------------   
*  (a) FD-TOA    Top-of-Atmosphere Radiative Fluxes
*  (b) FD-SRF    Surface Radiative Fluxes
*  (c) FD-PRF    Profiles  Radiative Flux Profiles (TOA and SRF, inclusive)
*  (d) FD-INP    Complete Input dataset which is only summarized in the 
*                other datasets
*  (e) FD-MPF    (Radiatively linearly averaged) Monthly-mean (hour) of FD-PRF
* -------------------------------------------------------------------------- 
*
*      All the above ISCCP-FD TOA, SRF, PRF, INP and MPF flux sub-datasets
* are written in big-endian, Integer*2, sequential-binary unformatted on
* UNIX machine with AIX OS. For those platforms (e.g., Intel PC) with 
* little-endian as default (Fortran) compiling option, one must explicitly 
* specify the big-endian option when doing compiling.
* 
*      The above regular global flux datasets, (a) to (d), have the same
* temporal resolution as ISCCP D1 (every three GMT's, 0, 3, .. 21). And (e) 
* is monthly mean of FD-PRF. They cover the time period from July 1983
* to the near-current. As of Septembet 2005, its ending time is extended to
* Dec. 2006 (from previous June, 2001) with 9701-0106 replaced by updated
* version because ISCCP-D1/D2 of 9801-0106 was updated to higher versions
* after last/first ISCCP-FD production in the end of 2002 (see ISCCP web
* site: "what's new?" at http://isccp.giss.nasa.gov/announcements.shtml)
* 
*      For a month (Integer*2 format), the total size of all the above
* datasets is about 0.6 GB, and separately,they have sizes of about
* 60, 70, 205,205 and 0.83 MB, respectively. For FD-MPF, the total dataset
* size for the currently available 18 years (216 months) is about 180 MB. 
* 
*     Appropriate scale factors for converting Integer*2 to physical Real*4 
* values are used in each read program. All the data file names for the above
* Integer*2 ISCCP-FD products, their Fortran read program names and the numbers
* of their output variables are as follows:  
*
* -------------------------------------------------------------------------
*    Secondary Integer*2 ISCCP-FD Product Names and their Read Program Names
* -------------------------------------------------------------------------- 
* Product Name | data file name  |    Name of       |  Number of output   
*                                |  Read Program    |     Parameters
* --------------------------------------------------------------------------  
* (a) FD-TOA  i2_toaii.YYMMDDHH  read_i2_FD_toa.f         20
* (b) FD-SRF  i2_srfii.YYMMDDHH  read_i2_FD_srf.f         25
* (c) FD-PRF  i2_prfii.YYMMDDHH  read_i2_FD_prf.f         81
* (d) FD-INP  i2_inpii.YYMMDDHH  read_i2_FD_inp.f    variable(max/av=281/128)
* (e) FD-MPF  i2_prfii.YYMM__..  read_i2_FD_mpf.f         81  (monthly mean)
* -------------------------------------------------------------------------- 
* 
* where "ii" = version. For 3-hourly data sets (a to d), YYMMDDHH =
* year-month-day-GMT, e.g., 99071521 is for year=1999, month=July,
* Day=15th and GMT=21, and for monthly mean (e), YYMM__.. = year-month,
* e.g., 9907__.. is for the monthly mean of July, 1999, where "__.." means 
* this monthly-mean product is derived by first averaging over a month for
* each GMT to eight monthly-hourly means and then averaging over these
* eight GMT's monthly-hourly means to monthly means (hour). 
* 
*     All the parameters from these integer*2 datasets may have slightly lower 
* precision than their original Real*4 [roughtly 0.5*1/(scale factor) for (a)
* to (d), and since (e) is averaged from (c), it may have approximately
* the same order of the precision as (c).
* 
* .......................................................................
*        WHEN USING THE ISCCP-FD DATA, PLEASE USE THE REFERENCE:
*
*	Zhang, y., w. b. Rossow, a. a. Lacis, v. Oinas, and
* m. i. Mishchenko (2004), Calculation of radiative fluxes from the
* surface to top of atmosphere based on ISCCP and other global data
* sets: Refinements of the radiative transfer model and the input data,
* j. Geophys. Res., 109, D19105, doi:10.1029/2003JD004457. 
***************************************************************************
*     The ISCCP FD-SRF RadFlux dataset gives the Upwelling and Downwelling,
* Shortwave (SW = 0.2-5.0 microns) and Longwave (LW = 5.0-200.0 microns)
* Radiative Fluxes at the surface, and a summary of the physical quantities
* used to calculate them.
*    All clear-sky fluxes are calculated under mean state for
* the atmosphere without clouds. The full-sky fluxes are weighted from
* fluxes calculated for ISCCP-defined 15 types of cloud condition and
* the clear-sky fluxes by their areas for each cell.  The "100% overcast"  
* fluxe exported from this program are reproduced from clear- and full-sky  
* fluxes using 
*       Overcast = [(full sky) - (Clear sky)*(1-CF)]/CF
* where CF is mean cloud cover, if CF > 0; If CF =0, i.e., the cell is   
* fully clear, all the three scenes have identical flux values since
* there is no cloud information available that can be used to calculate
* cloudy fluxes. Therefore, such overcast fluxes should be regarded as
* "apparent" to distinguish from the original overcast fluxes that
* have not been saved in order to save data storage space. 
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*    Read in global ISCCP FD-SRF (RadFlux surface flux) data and output 25 
* parameters in equal-area (EGA) or square (lat)-lon (equal-angle) maps. The 
* latter's longitudinal indices start from Dateline (SQD) or Greenwich line 
* (SQG).        
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
* Import:
*    Integer*2 collective global RadFlux surfae flux data file in 
* "i2_srfii.YYMMDDHH", where "ii" = version and YY/MM/DD/HH = 
* year/month/date/GMT, e.g.,99071521 is for year=1999, month=July, 
* Date=15th and GMT=21 
*.........................................................................
* Export:
*    25 parameters as listed below. Users have three options on output
* maps: EGA, SQD or SQG. The EGA has spatial resolution of equatorial
* 2.5 x 2.5 degree, about 280 km^2,written in real*4(6596), whose first 
* cell is for latitude of 90 s->87.5 s and longitude of 0 e->2.5 e. 
* Replicated from EGA map, both the SQD and SQG maps are of 2.5 x 2.5 degree,
* written in real*4(144,72). Their first latitudinal indices are for the cell
* of  90 s->87.5 s (South Pole), and their first longitudinal indices are for
* the cells of 180 w->177.5 w (Dateline) and 0 e->2.5 e(Greenwich line),
* respectively. Value -1000.00 is for no-data grid cells.
*     Output file names are in XXXXXXVV_YYMMDDHH, XXXXXXVV.YYMMDDHH,
* XXXXXXVV-YYMMDDHH for the above 3 different maps (EGA, SQD and SQG),
* repectively, where XXXXXX is output parameter name as listed and defined
* below, VV is version (='ii' now) and YYMMDDHH is from input file.
*......................................................................... 
*     The 25 parameters for each grid cell are:
*-----------------------------------------------------------------------
*Scale |Serial|     |                definition
*factor|number| name|
*-----------------------------------------------------------------------
*    10. (1) ps____  Surface pressure (mb);
*    10. (2) ts____  Surface skin temperature (Kelven);
* 10000. (3) al_srf  Surface broadband SW (0.2-5.0 micron) albedo (0-1);
* 10000. (4) em_srf  Surface broadband LW (5.0-200.0 micron) emissivity (0-1);
*  1000. (5) tlpwfl  Total column precipitable water for full sky (cm);
* 10000. (6) mu0___  Cosine solar zenith angle (0-1)
*    10. (7) tlo3__  Total column ozone (Dobson units);
*    10. (8) ta____  Surface air temperature (Kelven);
*  1000. (9) pws200  Precipitable water for 200-mb-thick layer covering 
*                    the surface for full sky (cm);
*  1000.(10) cf_m__  Mean cloud fraction (0-1);
*    10.(11) tau_m_  Mean cloud optical thickness;
*    10.(12) pc_m__  Mean cloud top  pressure (mb);
*    10.(13) pb_m__  Mean cloud base pressure (mb);
*    10.(14) sxdwbt  SW downwelling for full sky at surface (w/m^2);
*    10.(15) sxupbt  SW   upwelling for full sky at surface (w/m^2);
*    10.(16) txdwbt  LW downwelling for full sky at surface (w/m^2);
*    10.(17) txupbt  LW   upwelling for full sky at surface (w/m^2);
*    10.(18) srdbcr  SW downwelling for clear sky at surface (w/m^2);
*    10.(19) srubcr  SW   upwelling for clear sky at surface (w/m^2);
*    10.(20) trdbcr  LW downwelling for clear sky at surface (w/m^2);
*    10.(21) trubcr  LW   upwelling for clear sky at surface (w/m^2);
*---------------------------- Derived from above --------------------------
*       (22) sxdbcl  SW downwelling for 100% overcast sky at surface (w/m^2);
*       (23) sxubcl  SW   upwelling for 100% overcast sky at surface (w/m^2);
*       (24) txdbcl  LW downwelling for 100% overcast sky at surface (w/m^2);
*       (25) txubcl  LW   upwelling for 100% overcast sky at surface (w/m^2);
*..........................................................................
* Supplemental Subroutines: (1) eqar25(); (2) eq25ea()
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
* Users must specify (or comment out what they don't need):
*
* (1) uf_prd = the directory where input files are located
* (2) fn_prd = input file name="i2_srfii" as explained above
* (3) ufou   = the destination directory for output files 
*              (output parameter file name  = XXXXXXVV as explained above) 
* (4) iy/im/id/ih = year/month/date/GMT, e.g., 99/07/15/21, which defines
*                   YYMMDDHH for both the input and output files
* (5) iopout = choice for three different output maps:
*        0  output EGA map of Real*4(6596)   in XXXXXXVV_YYMMDDHH
*        1  output SQD map of Real*4(144,72) in XXXXXXVV.YYMMDDHH
*        2  output SQG map of Real*4(144,72) in XXXXXXVV-YYMMDDHH
* (6) prtlat/prtlon = latitude/longitude (prtlat=-90->90; prtlon=0->360 e)
*               of which the cell's 25 paramters can be printed out 
*               on screen for monitoring/checking
* (7) iwrite = 1, will write binary output files for a chosen Iopout
* (8) Lonlat = 1, will produce ceter longitude/latitude of all cells
*                 for a chosen map for users' convenience:
*                 elon(6596)/elat(6596) for equal-area map and
*                 alon(144,72)/alat(144,72) for equal-angle maps
**************************************************************************
%}
 %{ 
      parameter (iy=99, im=7, id=15, ih=21)
      parameter (iopout=0, iwrite=0, lonlat=0)
      parameter (undef=-1000.,iundef=-1000,
     &           prtlat=, prtlon=)
      parameter (Nsrf=21, nall=25,maxbox=6596, maxlat=72,maxlon=144) 
      dimension  srfflx(Nsrf,6596), eqrout(6596,nall),
     &           fluxrp(maxlon,maxlat), elon(maxbox),elat(maxbox),
     &           alon(144,72),alat(144,72), scfact(Nsrf)
      Integer*2  intsrf(Nsrf,6596)
      character  uf_prd*55,fn_prd*8,      fn*8,ft*8,
     &           ufou*55,  fn_srf(nall)*8
      logical    Lexist
      INTEGER*2 boxnmb,loneqr,lateqr,ncells  
      common/eqreqa/ boxnmb(144,72),loneqr(6596),lateqr(6596),     
     &               ncells(72),dlontb(72),xlatb(72),xlate(72)  
  %} 
%{
     sxdwbt_14=cell(1,8); 
     sxupbt_15=cell(1,8);
     txdwbt_16=cell(1,8);
     txupbt_17=cell(1,8);
     srdbcr_18=cell(1,8);
     srubcr_19=cell(1,8);
     trdbcr_20=cell(1,8);
     trubcr_21=cell(1,8);
     sxdbcl_22=cell(1,8);
     sxubcl_23=cell(1,8);
     txdbcl_24=cell(1,8);
     txubcl_25=cell(1,8);
%}
%{
     sxdwbt_14=zeros(1,8); 
     sxupbt_15=zeros(1,8);
     txdwbt_16=zeros(1,8);
     txupbt_17=zeros(1,8);
     srdbcr_18=zeros(1,8);
     srubcr_19=zeros(1,8);
     trdbcr_20=zeros(1,8);
     trubcr_21=zeros(1,8);
     sxdbcl_22=zeros(1,8);
     sxubcl_23=zeros(1,8);
     txdbcl_24=zeros(1,8);
     txubcl_25=zeros(1,8);
     %}
     

%{
     sxdwbt_14=zeros(8,31);
     sxupbt_15=zeros(8,31);
     txdwbt_16=zeros(8,31);
     txupbt_17=zeros(8,31);
     srdbcr_18=zeros(8,31);
     srubcr_19=zeros(8,31);
     trdbcr_20=zeros(8,31);
     trubcr_21=zeros(8,31);
     sxdbcl_22=zeros(8,31);
     sxubcl_23=zeros(8,31);
     txdbcl_24=zeros(8,31);
     txubcl_25=zeros(8,31);
%}
%{
     sxdwbt_14=zeros(8,31,12);
     sxupbt_15=zeros(8,31,12);
     txdwbt_16=zeros(8,31,12);
     txupbt_17=zeros(8,31,12);
     srdbcr_18=zeros(8,31,12);
     srubcr_19=zeros(8,31,12);
     trdbcr_20=zeros(8,31,12);
     trubcr_21=zeros(8,31,12);
     sxdbcl_22=zeros(8,31,12);
     sxubcl_23=zeros(8,31,12);
     txdbcl_24=zeros(8,31,12);
     txubcl_25=zeros(8,31,12);



     sxdwbt_14=zeros(144,72,8,31,12);
     sxupbt_15=zeros(144,72,8,31,12);
     txdwbt_16=zeros(144,72,8,31,12);
     txupbt_17=zeros(144,72,8,31,12);
     srdbcr_18=zeros(144,72,8,31,12);
     srubcr_19=zeros(144,72,8,31,12);
     trdbcr_20=zeros(144,72,8,31,12);
     trubcr_21=zeros(144,72,8,31,12);

%}

%{

     sxdbcl_22=zeros(8,31,12,27);
     sxubcl_23=zeros(8,31,12,27);
     txdbcl_24=zeros(8,31,12,27);
     txubcl_25=zeros(8,31,12,27);
%}
   
     ssr_=zeros(144,72,26);
     ssrc_=zeros(144,72,26);
     ssrd_=zeros(144,72,26);
     str_=zeros(144,72,26);
     strc_=zeros(144,72,26);
     strd_=zeros(144,72,26);

tic;
%{
for iy=00:09
for im=1:12
for id=1:31
%}
for iy=1984:2009
    
    
    
     sxdwbt_14=zeros(144,72,8,31,12);
     sxupbt_15=zeros(144,72,8,31,12);
     txdwbt_16=zeros(144,72,8,31,12);
     txupbt_17=zeros(144,72,8,31,12);
     srdbcr_18=zeros(144,72,8,31,12);
     srubcr_19=zeros(144,72,8,31,12);
     trdbcr_20=zeros(144,72,8,31,12);
     trubcr_21=zeros(144,72,8,31,12);
    
    
for im=1:12
for id=1:31
for ih=0:3:21
%iy=01; im=01; %id=01; %ih=00;


     
     
     
     
     
     
iopout=2; iwrite=0; lonlat=1;
undef=-1000; iundef=-1000; prtlon=285.0; prtlat=43.0;
nsrf=21; nall=25; maxbox=6596;  maxlat=72; maxlon=144;
fluxrp=ones(maxlon,maxlat);
%srfflx=ones(nsrf,maxbox); eqrout=ones(maxbox,nall);
%elon=ones(1,maxbox); elat=ones(1,maxbox);
%alon=ones(144,72); alat=ones(144,72);

 boxnmb=ones(144,72); loneqr=ones(1,6596); lateqr=ones(1,6596);
 ncells=ones(1,72); dlontb=ones(1,72); xlatb=ones(1,72); xlate=ones(1,72);
 
 

 %{
 switch fix(iy/10)
    case 8
        iy_=num2str(1900+iy);
    case 9
        iy_=num2str(1900+iy);
    case 0
        iy_=num2str(2000+iy);
end

 
switch im
    case 1
        im_='January';
    case 2
        im_='February';
    case 3
        im_='March';
    case 4
        im_='April';
    case 5
        im_='May';
    case 6
        im_='June';
    case 7
        im_='July';
    case 8
        im_='August';
    case 9
        im_='September';
    case 10
        im_='October';
    case 11
        im_='November';
    case 12
        im_='December';
end

switch mod(id,10)
    case 1
        id_=strcat(num2str(id),'st');
    case 2
        id_=strcat(num2str(id),'nd');
    case 3
        id_=strcat(num2str(id),'rd');
    case 4
        id_=strcat(num2str(id),'th');
    case 5
        id_=strcat(num2str(id),'th');
    case 6
        id_=strcat(num2str(id),'th');
    case 7
        id_=strcat(num2str(id),'th');
    case 8
        id_=strcat(num2str(id),'th');
    case 9
        id_=strcat(num2str(id),'th');
    case 0
        id_=strcat(num2str(id),'th');
end
%} 
 
%c...................................................................
%{
      data scfact/2*10., 2*10000.; 1000.;10000.; 2*10.; 2*1000.; 11*10./
%c:    input file name:
      data fn_prd/'i2_srfii'/
%c:    input directory:
      data uf_prd/'/aos/home/syin/Desktop/        .          '/,
%c                            11  15                                  
     &     ipprd1/10/, ipprd2/15/  
%c:    output destination:
      data ufou /'/aos/home/syin/Desktop/out_/        .          '/,
%c                                   19   24  
     &     ipou1/19/, ipou2/24/
%c:    25 output parameters' file names -- users may freely change them:
      data  fn_srf/
     1'ps____  ','ts____  ','al_srf  ','em_srf  ','tlpwfl  ','mu0___  ',
     2'tlo3__  ','ta____  ','pws200  ','cf_m__  ','tau_m_  ','pc_m__  ',
     3'pb_m__  ','sxdwbt  ','sxupbt  ','txdwbt  ','txupbt  ','srdbcr  ',
     4'srubcr  ','trdbcr  ','trubcr  ','sxdbcl  ','sxubcl  ','txdbcl  ',
     5'txubcl  '/
%}     
scfact=[10,10,10000,10000,1000,10000,10,10,1000,1000,10,10,10,10,10,10,10,10,10,10,10];
fn_prd='i2_srfii';
%uf_prd='/aos/home/syin/Desktop/.';
%ipprd1=10;
%ipprd2=15;
%ufou='/aos/home/syin/Desktop/out_/.';
%ipou1=19;
%ipou2=24;
fn_srf={'ps____  ','ts____  ','al_srf  ','em_srf  ','tlpwfl  ','mu0___  ','tlo3__  ','ta____  ','pws200  ','cf_m__  ','tau_m_  ','pc_m__  ','pb_m__  ','sxdwbt  ','sxupbt  ','txdwbt  ','txupbt  ','srdbcr  ','srubcr  ','trdbcr  ','trubcr  ','sxdbcl  ','sxubcl  ','txdbcl  ','txubcl  '};
  
%c--------------------------------------------------------------------
    %CALL eqar25(2.5,2.5)
    [dlontb,ncells,loneqr,lateqr,xlatb,xlate,boxnmb]=eqar25(2.5,2.5);
%c--
%c   (1) Read in the collective global SRF RadFlux data file i2_srfii.YYMMDDHH
%c--
%c:    define fn='i2_srfii' and ft='YYMMDDHH':
      FN=fn_prd;
%      write(ft(1:8),'(4i2.2)')iy,im,id,ih
       fn='i2_srfii';
       ft=strcat(num2str(mod(iy,100),'%2.2i'),num2str(im,'%2.2i'),num2str(id,'%2.2i'),num2str(ih,'%2.2i'));
%c:    initialization:
    srfflx=ones(nsrf,maxbox);
    intsrf=ones(nsrf,maxbox);
      for ibx=1:maxbox
      for ip =1:nsrf
         srfflx(ip,ibx)=undef;
         intsrf(ip,ibx)=iundef;
      end 
      end 
%c:       define and read in input file: 
%cc       uf_prd(ipprd1+1:ipprd1+4)=ft(1:4)
%         uf_prd(ipprd2+1:ipprd2+19)=fn(1:8)//'.'//ft(1:8);
          uf_prd(1:17)=[fn(1:8),'.',ft(1:8)];
%         inquire(file=uf_prd, exist=Lexist);
%          file = dir('/aos/home/syin/MATLAB workspace/project3/ISCCP/Main and GD/*');
%          file = dir('/aos/home/syin/MATLAB workspace/project3/ISCCP/i2_srfii.2001/*');          
%         file = dir('/aos/home/syin/MATLAB workspace/project3/ISCCP/i2_srfii.from2000to2009/*');
         file = dir('/aos/home/syin/MATLAB workspace/project4/Bias between ERA Interim and ISCCP/i2_srfii.from1983to2009/*');
%         file = dir('/aos/home/syin/Downloads/*');
      x=ones(1,numel(file));
    for i=1:numel(file)
      x(i)=strcmp(uf_prd,file(i).name);
    end 
       if sum(x)==0
%       if Lexist~=true
        %   print( *,'................................................')
        %   print( *,'----   NON-EXISTENT FILE: ',uf_prd )
            fprintf(1,'................................................\n');
            fprintf(1,['----   NON-EXISTENT FILE: ',uf_prd,'\n']);
%            break;
       else
        %  print( *,'++++ 11 ++++ UF_PRD = ',UF_PRD)
           fprintf(1,['++++ 11 ++++ UF_PRD = ',uf_prd,'\n']);
        %  OPEN(UNIT=11,FILE=Uf_prd,
           fid_11=fopen(uf_prd,'r','b');
   %  *  STATUS='old',ACCESS='sequential',FORM='UNFORMATTED',
   %  *  IOSTAT=IOS)
   %       read(11)intsrf
        intsrf=ones(nsrf,maxbox);
          intsrf_=fread(fid_11,'int16');  
    %     close (11)
          status=fclose(fid_11);
       end
%c
    m=3;
    for j=1:maxbox
        for i=1:nsrf
           intsrf(i,j)=intsrf_(m);
           m=m+1;
        end
    end
      for  i = 1:maxbox                      %DO 222
      for  k = 1:nsrf                        %DO 222
         if intsrf(k,i)>iundef
            srfflx(k,i)=intsrf(k,i)/scfact(k);
         else
            srfflx(k,i)=-1000;
         end 
                                %222 CONTINUE
      end
      end
%c--
%c   (2) Get equal-area maps for all 25 output parameters
%c--
%c:    Directly get first 21 parameters
  eqrout=ones(maxbox,nall);
      for  ip = 1:nall                           %DO 121
         for ibx = 1:maxbox
            if ip<=nsrf
               eqrout(ibx,ip)=srfflx(ip,ibx); 
            else
               eqrout(ibx,ip)=undef;
            end
         end 
                                     %121 CONTINUE
      end
%c: Calculate apparent overcast fluxes when a cell is not fully clear
%c and equate it to full-/clear-sky when the cell is fully clear sky (in this
%c case, no cloudy flux can be produced since no cloud info is available)
%{
       for  ibx = 1:maxbox                    %DO 131
            CF=eqrout(ibx,10);
         if CF>0
           if((srfflx(14,ibx)>undef)&&(srfflx(18,ibx)>undef))
           eqrout(ibx,22)=(srfflx(14,ibx)-srfflx(18,ibx)*(1-CF))/CF;
           end
           if((srfflx(15,ibx)>undef)&&(srfflx(19,ibx)>undef))
           eqrout(ibx,23)=(srfflx(15,ibx)-srfflx(19,ibx)*(1-CF))/CF;
           end
           if((srfflx(16,ibx)>undef)&&(srfflx(20,ibx)>undef))
           eqrout(ibx,24)=(srfflx(16,ibx)-srfflx(20,ibx)*(1-CF))/CF;
           end
           if((srfflx(17,ibx)>undef)&&(srfflx(21,ibx)>undef))
           eqrout(ibx,25)=(srfflx(17,ibx)-srfflx(21,ibx)*(1-CF))/CF;
           end
         else if CF==0
           if(srfflx(14,ibx)>undef)
           eqrout(ibx,22)=srfflx(14,ibx);
           end
           if(srfflx(15,ibx)>undef)
           eqrout(ibx,23)=srfflx(15,ibx);
           end
           if(srfflx(16,ibx)>undef)
           eqrout(ibx,24)=srfflx(16,ibx);
           end
           if(srfflx(17,ibx)>undef)
           eqrout(ibx,25)=srfflx(17,ibx);
           end 
             end
         end
                                 %131 CONTINUE
      end
%}
%c--
%c   (3) output for iopout = 0, or 1, or 2
%c--
   elat=ones(1,maxbox);
   elon=ones(1,maxbox);
      if(iopout==0)   %! output EGA 
        if(lonlat==1)
%c:      produce center longitude/latitude of each cell for EGA map
           for ibx=1:maxbox
                lat=lateqr(ibx);
              elat(ibx)=((lat)-0.5)*2.5 - 90;
              elon(ibx)=(loneqr(ibx)-0.5)*dlontb(lat);
%cc            if(loneqr(ibx).eq.ncells(lat))
%cc   &     print  *,'ibx/elat/elon=',ibx,'/',elat(ibx),'/',elon(ibx)
           end
        end
%c
         lat=floor((prtlat+90)/2.5) + 1;
         if lat>maxlat
             lat=maxlat;
         end
         lon=floor(prtlon/dlontb(lat))+1;
         if lon>ncells(lat)
             lon=ncells(lat);
         end
         jbx=boxnmb((lon),lat);
        for  ip = 1:nall                         %DO 221
%c:       printout parameters of the cell for the designated longitude/latitude
%         write(6,6221)fn_srf(ip)(1:6),eqrout(jbx,ip),prtlon,prtlat,jbx
% 6221    format(a6,' = ',f10.4,' for prtlon/prtlat = ',f6.2,'/',f6.2,
%     &      ' & eq-area box #=',i4)
%          fwrite(1,[fn_srf{ip}(1:6),eqrout(jbx,ip),prtlon,prtlat,jbx]);  
         fprintf(1,'%6s = %10.4f for prtlon/prtlat = %6.2f/%6.2f & eq-area box #= %4i \n',fn_srf{ip}(1:6),eqrout(jbx,ip),prtlon,prtlat,jbx);
          if(iwrite==1)
%c:          write out binary files for 25 paramters
            fn=[fn_srf{ip}(1:6),fn_prd(7:8)];
%cc          if(iyymmo.eq.1)ufou(ipou1+1:ipou1+4)=ft(1:4)
            ufou(1:21)=[fn(1:8),'_',ft(1:8),'.txt'];
   %         print( *,'++++ 21 ++++ ufou = ',ufou);
             fprintf(1,['++++ 21 ++++ ufou = ',ufou]);
   %         OPEN(UNIT=21,FILE=ufou,
            fid_21=fopen(ufou,'a');
   %  *      STATUS='unknown',ACCESS='sequential',FORM='UNFORMATTED',
   %  *      IOSTAT=IOS)
   %         WRITE(21)(eqrout(ibx,ip),ibx=1,maxbox)
             for ibx=1:maxbox
             fprintf(fid_21,'%d',eqrout(ibx,ip));
             end
   %         close (21)
            status=fclose(fid_21);
          end 
     continue;       %221 CONTINUE
        end
   alat=ones(maxlon,maxlat);
   alon=ones(maxlon,maxlat);
      else if iopout==1  %! output output SQD
%c:      produce center longitude/latitude of each cell for SQD 
        if(lonlat==1)
           for j = 1:maxlat
           for i = 1:maxlon
              alat(i,j)=(j-0.5)*2.5 - 90;
              alon(i,j)=(i-0.5)*2.5 +180;
              if alon(i,j)>360
                  alon(i,j)=alon(i,j)-360;
%cc            if(i.eq.maxlon)
%cc   &     print *,'j,i,alat,alon=',j,'/',i,'/',alat(i,j),'/',alon(i,j)
              end 
           end 
           end 
        end
         lat=(prtlat+90)/2.5 + 1;
         if lat>maxlat
             lat=maxlat;
         end
         lon=prtlon/2.5+1;
         if(lon>144)
             lon=144;
         end
         lon=lon+72;
         if(lon>144)
             lon=lon-144;
         end
        for  ip = 1:nall                 %DO 231
    %CALL eq25ea(1,eqrout(1,ip),fluxrp)
     [boxnmb,eqangl]=eq25ea(1,eqrout(1,ip),fluxrp);
%c:       printout parameters of the cell for the designated longitude/latitude
%         write(6,6231)fn_srf(ip)(1:6),fluxrp(lon,(lat)),
%     &                prtlon,prtlat,lon,(lat)
% 6231    format(a6,'=',f10.4,' for prtlon/prtlat = ',f6.2,'/',f6.2,
%     &      ' & SQD index lon/lat= ',i3,'/',i2)
         
           if(iwrite==1)
%c:          write out binary files for 25 paramters
            fn=[fn_srf{ip}(1:6),fn_prd(7:8)];
%cc          if(iyymmo.eq.1)ufou(ipou1+1:ipou1+4)=ft(1:4)
            ufou(1:21)=[fn(1:8),'.',ft(1:8),'.txt'];
    %        print( *,'++++ 21 ++++ ufou = ',ufou);
             fprint(1,'++++ 21 ++++ ufou = ',ufou);
    %        OPEN(UNIT=21,FILE=ufou,
             fid_21=fopen(ufou,'a');
    % *      STATUS='unknown',ACCESS='sequential',FORM='UNFORMATTED',
    % *      IOSTAT=IOS)
    %        WRITE(21)fluxrp
             fprintf(fid_21,'%d',fluxrp);
    %       close (21)
            status=fclose(fid_21);
          end 
     continue;                           %231 CONTINUE
        end
         else if(iopout==2)  %! output SQG
%c:      produce center longitude/latitude of each cell for SQG 
        if(lonlat==1)
           for j = 1:maxlat
           for i = 1:maxlon
              alat(i,j)=(j-0.5)*2.5 - 90;
              alon(i,j)=(i-0.5)*2.5; 
%cc            if(i.eq.maxlon)
%cc   &     print *,'j,i,alat,alon=',j,'/',i,'/',alat(i,j),'/',alon(i,j)
           end 
           end 
        end 
         lat=floor((prtlat+90.)/2.5) + 1;
         if lat>maxlat
             lat=maxlat;
         end
         lon=floor(prtlon/2.5)+1;
         if lon>144
             lon=144;
         end
        for  ip = 14:21                  %DO 241
     %CALL eq25ea(0,eqrout(1,ip),fluxrp)
     [boxnmb,fluxrp]=eq25ea(0,eqrout(:,ip),fluxrp);
     
     
    lat_=-88.75:2.5:88.75; 
    
     switch ip 
         
         case 14
             sxdwbt_14(:,:,ih/3+1,id,im)=fluxrp;
%             sxdwbt_14(ih/3+1,id,im,iy-1983)=area_mean(squeeze(mean(fluxrp)),lat_);
%             sxdwbt_14(ih/3+1,id,im,iy)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            sxdwbt_14(ih/3+1,id,im)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            sxdwbt_14(ih/3+1,id)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            sxdwbt_14(ih/3+1)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            sxdwbt_14(ih/3+1)=squeeze(mean(mean(fluxrp)));    
%             sxdwbt14GD(fluxrp,iy_,im_,id_,ih);
%             sxdwbt_14=fluxrp;
         case 15
             sxupbt_15(:,:,ih/3+1,id,im)=fluxrp;
%             sxupbt_15(ih/3+1,id,im,iy-1983)=area_mean(squeeze(mean(fluxrp)),lat_);
%             sxupbt_15(ih/3+1,id,im,iy)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            sxupbt_15(ih/3+1,id,im)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            sxupbt_15(ih/3+1,id)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            sxupbt_15(ih/3+1)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%             sxupbt_15(ih/3+1)=squeeze(mean(mean(fluxrp)));
%            sxupbt15GD(fluxrp,iy_,im_,id_,ih);
%             sxupbt_15=fluxrp;
         case 16
             txdwbt_16(:,:,ih/3+1,id,im)=fluxrp;
%             txdwbt_16(ih/3+1,id,im,iy-1983)=area_mean(squeeze(mean(fluxrp)),lat_);
%             txdwbt_16(ih/3+1,id,im,iy)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            txdwbt_16(ih/3+1,id,im)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            txdwbt_16(ih/3+1,id)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            txdwbt_16(ih/3+1)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%             txdwbt_16(ih/3+1)=squeeze(mean(mean(fluxrp)));
%             txdwbt16GD(fluxrp,iy_,im_,id_,ih);
%             txdwbt_16=fluxrp;
         case 17
             txupbt_17(:,:,ih/3+1,id,im)=fluxrp;
%             txupbt_17(ih/3+1,id,im,iy-1983)=area_mean(squeeze(mean(fluxrp)),lat_);
%             txupbt_17(ih/3+1,id,im,iy)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            txupbt_17(ih/3+1,id,im)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            txupbt_17(ih/3+1,id)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            txupbt_17(ih/3+1)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%             txupbt_17(ih/3+1)=squeeze(mean(mean(fluxrp)));
%             txupbt17GD(fluxrp,iy_,im_,id_,ih);
%             txupbt_17=fluxrp;
         case 18
             srdbcr_18(:,:,ih/3+1,id,im)=fluxrp;
%             srdbcr_18(ih/3+1,id,im,iy-1983)=area_mean(squeeze(mean(fluxrp)),lat_);
%             srdbcr_18(ih/3+1,id,im,iy)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            srdbcr_18(ih/3+1,id,im)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            srdbcr_18(ih/3+1,id)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            srdbcr_18(ih/3+1)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%             srdbcr_18(ih/3+1)=squeeze(mean(mean(fluxrp)));
%             srdbcr18GD(fluxrp,iy_,im_,id_,ih);
%             srdbcr_18=fluxrp;
         case 19
             srubcr_19(:,:,ih/3+1,id,im)=fluxrp;
%             srubcr_19(ih/3+1,id,im,iy-1983)=area_mean(squeeze(mean(fluxrp)),lat_);
%             srubcr_19(ih/3+1,id,im,iy)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            srubcr_19(ih/3+1,id,im)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            srubcr_19(ih/3+1,id)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            srubcr_19(ih/3+1)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%             srubcr_19(ih/3+1)=squeeze(mean(mean(fluxrp)));
%             srubcr19GD(fluxrp,iy_,im_,id_,ih);       
%             srubcr_19=fluxrp;
         case 20
             trdbcr_20(:,:,ih/3+1,id,im)=fluxrp;
%             trdbcr_20(ih/3+1,id,im,iy-1983)=area_mean(squeeze(mean(fluxrp)),lat_);
%             trdbcr_20(ih/3+1,id,im,iy)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            trdbcr_20(ih/3+1,id,im)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            trdbcr_20(ih/3+1,id)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            trdbcr_20(ih/3+1)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%             trdbcr_20(ih/3+1)=squeeze(mean(mean(fluxrp)));
%             trdbcr20GD(fluxrp,iy_,im_,id_,ih);    
%             trdbcr_20=fluxrp;
         case 21
             trubcr_21(:,:,ih/3+1,id,im)=fluxrp;
%             trubcr_21(ih/3+1,id,im,iy-1983)=area_mean(squeeze(mean(fluxrp)),lat_);
%             trubcr_21(ih/3+1,id,im,iy)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            trubcr_21(ih/3+1,id,im)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            trubcr_21(ih/3+1,id)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            trubcr_21(ih/3+1)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%             trubcr_21(ih/3+1)=squeeze(mean(mean(fluxrp)));
%             trubcr21GD(fluxrp,iy_,im_,id_,ih);   
%             trubcr_21=fluxrp;
%         case 22
%             sxdbcl_22(ih/3+1,id,im,iy)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            sxdbcl_22(ih/3+1,id,im)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            sxdbcl_22(ih/3+1,id)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            sxdbcl_22(ih/3+1)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            sxdbcl_22(ih/3+1)=fluxrp;
%            sxdbcl22GD(fluxrp,iy_,im_,id,ih);             
%         case 23
%             sxubcl_23(ih/3+1,id,im,iy)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            sxubcl_23(ih/3+1,id,im)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            sxubcl_23(ih/3+1,id)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            sxubcl_23(ih/3+1)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            sxubcl_23(ih/3+1)=fluxrp;
%            sxubcl23GD(fluxrp,iy_,im_,id,ih);             
%         case 24
%             txdbcl_24(ih/3+1,id,im,iy)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            txdbcl_24(ih/3+1,id,im)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            txdbcl_24(ih/3+1,id)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            txdbcl_24(ih/3+1)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            txdbcl_24(ih/3+1)=fluxrp;
%            txdbcl24GD(fluxrp,iy_,im_,id,ih);             
%         case 25
%             txubcl_25(ih/3+1,id,im,iy)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            txubcl_25(ih/3+1,id,im)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            txubcl_25(ih/3+1,id)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            txubcl_25(ih/3+1)=0.25*(fluxrp(lon,lat)+fluxrp(lon,lat-1)+fluxrp(lon-1,lat)+fluxrp(lon-1,lat-1));
%            txubcl_25(ih/3+1)=fluxrp;
%            txubcl25GD(fluxrp,iy_,im_,id,ih);                 
     
     
%     if ip==14,%15,%16,%17,%18,%19,%20,%21,%22,%23,%24,%25 
%     sxdwbt14GD(fluxrp,iy,im,id,ih);     
%     sxupbt15GD(fluxrp,iy,im,id,ih);
%     txdwbt16GD(fluxrp,iy,im,id,ih);
%     txupbt17GD(fluxrp,iy,im,id,ih);
%     srdbcr18GD(fluxrp,iy,im,id,ih);
%     srubcr19GD(fluxrp,iy,im,id,ih);
%     trdbcr20GD(fluxrp,iy,im,id,ih);
%     trubcr21GD(fluxrp,iy,im,id,ih);
%     sxdbcl22GD(fluxrp,iy,im,id,ih);
%     sxubcl23GD(fluxrp,iy,im,id,ih);
%     txdbcl24GD(fluxrp,iy,im,id,ih);
%     txubcl25GD(fluxrp,iy,im,id,ih);     
     end
     
%c:       printout parameters of the cell for the designated longitude/latitude
%        write(6,6241)fn_srf(ip)(1:6),fluxrp(lon,(lat)),
%     &               prtlon,prtlat,lon,(lat)
% 6241   format(a6,'=',f10.4,' for prtlon/prtlat = ',f6.2,'/',f6.2,
%     &      ' & SQG index lon/lat= ',i3,'/',i2)


%         fprintf(1,'%6s = %10.4f for prtlon/prtlat = %6.2f/ %6.2f & sqg index lon/lat= %3i/ %2i \n',fn_srf{ip}(1:6),fluxrp(lon,(lat)),prtlon,prtlat,lon,(lat));   

         
         if(iwrite==1)
%c:          write out binary files for 25 paramters
            fn=[fn_srf{ip}(1:6),fn_prd(7:8)];
%cc          if(iyymmo.eq.1)ufou(ipou1+1:ipou1+4)=ft(1:4)
   %         ufou(ipou2+1:ipou2+19)=fn(1:8)//'-'//ft(1:8);
             ufou(1:21)=[fn(1:8),'-',ft(1:8),'.txt'];
   %         print( *,'++++ 41 ++++ ufou = ',ufou);
             fprintf(1,'++++ 41 ++++ ufou = ',ufou);
   %         OPEN(UNIT=41,FILE=ufou,
             fid_41=fopen(ufou,'a');
   %  *      STATUS='unknown',ACCESS='sequential',FORM='UNFORMATTED',
   %  *      IOSTAT=IOS)
   %         WRITE(41)fluxrp
             fprintf(fid_41,'%d',fluxrp);
   %        close (41)
            status=fclose(fid_41);
         end 
                                   %241 CONTINUE
        end 
             end
          end
      end

      
end
end
end

for i=1:144
    for j=1:72
        
ssr=squeeze(nanmean((sxdwbt_14(i,j,:,:,:)-sxupbt_15(i,j,:,:,:)),3));
ssr_(i,j,iy-1983)=squeeze(nanmean(ssr(ssr~=0)));
ssrc=squeeze(nanmean((srdbcr_18(i,j,:,:,:)-srubcr_19(i,j,:,:,:)),3));
ssrc_(i,j,iy-1983)=squeeze(nanmean(ssrc(ssrc~=0)));
ssrd=squeeze(nanmean(sxdwbt_14(i,j,:,:,:),3));
ssrd_(i,j,iy-1983)=squeeze(nanmean(ssrd(ssrd~=0)));
str=squeeze(nanmean((txdwbt_16(i,j,:,:,:)-txupbt_17(i,j,:,:,:)),3));
str_(i,j,iy-1983)=squeeze(nanmean(str(str~=0)));
strc=squeeze(nanmean((trdbcr_20(i,j,:,:,:)-trubcr_21(i,j,:,:,:)),3));
strc_(i,j,iy-1983)=squeeze(nanmean(strc(strc~=0)));
strd=squeeze(nanmean(txdwbt_16(i,j,:,:,:),3));
strd_(i,j,iy-1983)=squeeze(nanmean(strd(strd~=0)));
    end
end




end
%{
end
end
end
%}     

%{      
ssr14_15GD(sxdwbt_14,sxupbt_15,iy_,im_,id_,ih);
str16_17GD(txdwbt_16,txupbt_17,iy_,im_,id_,ih);
ssrc18_19GD(srdbcr_18,srubcr_19,iy_,im_,id_,ih);
strc20_21GD(trdbcr_20,trubcr_21,iy_,im_,id_,ih);
%}
      



  save ('/aos/home/syin/MATLAB workspace/project4/workspace.mat');
%{
switch rem(iy,10)
    case 8
        iy_=num2str(1900+iy);
    case 9
        iy_=num2str(1900+iy);
    case 0
        iy_=num2str(2000+iy);
end

switch im
    case 1
        im_='January';
    case 2
        im_='February';
    case 3
        im_='March';
    case 4
        im_='April';
    case 5
        im_='May';
    case 6
        im_='June';
    case 7
        im_='July';
    case 8
        im_='August';
    case 9
        im_='September';
    case 10
        im_='October';
    case 11
        im_='November';
    case 12
        im_='December';
end
%}    

%{
sxdwbt14AM(sxdwbt_14);
sxupbt15AM(sxupbt_15);
txdwbt16AM(txdwbt_16);
txupbt17AM(txupbt_17);
srdbcr18AM(srdbcr_18);
srubcr19AM(srubcr_19);
trdbcr20AM(trdbcr_20);
trubcr21AM(trubcr_21);
sxdbcl22AM(sxdbcl_22);
sxubcl23AM(sxubcl_23);
txdbcl24AM(txdbcl_24);
txubcl25AM(txubcl_25);
%}
%{
sxdwbt14MM(sxdwbt_14);
sxupbt15MM(sxupbt_15);
txdwbt16MM(txdwbt_16);
txupbt17MM(txupbt_17);
srdbcr18MM(srdbcr_18);
srubcr19MM(srubcr_19);
trdbcr20MM(trdbcr_20);
trubcr21MM(trubcr_21);
sxdbcl22MM(sxdbcl_22);
sxubcl23MM(sxubcl_23);
txdbcl24MM(txdbcl_24);
txubcl25MM(txubcl_25);
%}

%{

sxdwbt14DM(sxdwbt_14,iy_,im_);
sxupbt15DM(sxupbt_15,iy_,im_);
txdwbt16DM(txdwbt_16,iy_,im_);
txupbt17DM(txupbt_17,iy_,im_);
srdbcr18DM(srdbcr_18,iy_,im_);
srubcr19DM(srubcr_19,iy_,im_);
trdbcr20DM(trdbcr_20,iy_,im_);
trubcr21DM(trubcr_21,iy_,im_);
sxdbcl22DM(sxdbcl_22,iy_,im_);
sxubcl23DM(sxubcl_23,iy_,im_);
txdbcl24DM(txdbcl_24,iy_,im_);
txubcl25DM(txubcl_25,iy_,im_);


%}


%{
ssr14_15DC(sxdwbt_14,sxupbt_15,iy_,im_,id_);
str16_17DC(txdwbt_16,txupbt_17,iy_,im_,id_);
ssrc18_19DC(srdbcr_18,srubcr_19,iy_,im_,id_);
strc20_21DC(trdbcr_20,trubcr_21,iy_,im_,id_);
sxdwbt14DC(sxdwbt_14,iy_,im_,id_);
txdwbt16DC(txdwbt_16,iy_,im_,id_);
%}


%{
sxupbt15DC(sxupbt_15,iy,im,id);
txupbt17DC(txupbt_17,iy,im,id);
srdbcr18DC(srdbcr_18,iy,im,id);
srubcr19DC(srubcr_19,iy,im,id);
trdbcr20DC(trdbcr_20,iy,im,id);
trubcr21DC(trubcr_21,iy,im,id);
sxdbcl22DC(sxdbcl_22,iy,im,id);
sxubcl23DC(sxubcl_23,iy,im,id);
txdbcl24DC(txdbcl_24,iy,im,id);
txubcl25DC(txubcl_25,iy,im,id);
%}


toc;
%{
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
      Function eqar25(dlon,dlat)
      dlon=0.5;dlat=0.5;  
      lonmin=1;latmin=1;
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
      lonmax=360.0/dlon;                                                 
      latmax=180.0/dlat;                                                 
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
     fprintf(1,'  dlat   latmax  area eq cell %8.2f %5i %8.1f',[dlat,latmax,aecell]);
%c         
   xlatb=ones(1,latmax);   
   xlate=ones(1,latmax);
   ncells=ones(1,latmax);
      for  (lat)=latmin:latmax              %DO 100                           
           xlatb((lat))=dlat*((lat)-1)-90.0;
           xlate((lat))=xlatb((lat))+dlat;                              
           if xlate((lat))>90
               xlate((lat))=90;
           end
           if xlatb((lat))>90
               xlatb((lat))=n:latmax;              %DO 100                           
               xlatb((lat))=dlat*((lat)-1)-90.0;
           xlate((lat))=xlatb((lat))+dlat;
           end
           if xlate((lat))>90
               xlate((lat))=90;
           end
           if xlatb((lat))>90
               xlatb((lat))=90;         
           end
      rlatb=twopir*xlatb((lat));                     
      rlate=twopir*xlate((lat));                                        
      htb=re*sin(rlatb);                                        
      hte=re*sin(rlate);                                                
      htzone=ABS(hte-htb);                                               
      azone=twopi*re*htzone;                                             
      cells=azone/aecell;                                             
      ncells((lat))=cells+0.5;
   continue;                              %100 CONTINUE  
      end
%c                                                                       
      latlft=latmin;                                          
      if  mod(latmax-latmin+1,2)==1 
           latrit=latmax-1;                                            
      else                                                              
           latrit=latmax;
      end                                                             
      for  (lat)=latmin:latmax         %DO 200                                  
           if ( latlft>=latrit )     %GOTO 201  
               continue;             %201 CONTINUE
           else
           ncells( latrit )=ncells( latlft );                            
           xlatb( latrit )=-xlate( latlft );                            
           xlate( latrit )=-xlatb( latlft );
           latlft=latlft+1;                                              
           latrit=latrit-1;     
           end
   continue;                        %200 CONTINUE
      end                                                        
%c                                                                       
%c     WRITE(31,3000)                                                    
% 3000 FORMAT(1X,'(lat)',5X,'xlatb',5X,'xlate',5X,'ncells',5X,'dlontb')   
%  dlontb=ones(latmax);
      fprintf(fid_31,'(lat)     xlatb     xlate     ncells     dlontb',[]);
      for  (lat)=latmin:latmax        %DO 300                                 
      if  ncells((lat))>0                                  
           dlontb((lat))=360./ncells((lat));
           else                                                              
           dlontb((lat))=360.0;                                          
      end                                                             
      if xlate((lat))>90
          xlate((lat))=90;
%c     WRITE(21,301) (lat),xlatb((lat)),xlate((lat)),                      
%c    &          ncells((lat)),dlontb((lat))                             
%  301 FORMAT(1X,'(lat)',I4,1X,'XLAT',2F8.1,2X,'CELLS',                    
%     1  I8,2X,'DLON',F8.1)    
      fprintf(fid_21,'(lat) xlat  cells dlon %4i %8.1f %8.1f %8i %8.1f',[(lat),xlatb((lat)),xlate((lat)),ncells(lat),dlontb(lat)]);
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
      function [ioptn,arybox,eqangl]= eq25ea(ioptn,arybox,eqangl)
      arybox=ones(1,6596);eqangl=ones(144,72);
      boxnmb=1; loneqr=1; lateqr=1; ncells=1;     
      boxnmb=ones(144,72); loneqr=ones(1,6596); lateqr=ones(1,6596);     
      ncells=ones(1,72); dlontb=ones(1,72); xlatb=ones(1,72); xlate=ones(1,72);
      for  j=1:72           %DO 101                                          
      for  i=1:144          %DO 101                                          
          eqangl(i,j)=-1000;                                            
      continue;                %101 CONTINUE    
      end
      end
      for  j=1:72             %DO 111                                        
          (lat)=j;                                                        
      for  i=1:144            %DO 111                                        
          if(ioptn==0)                                           
             iadj=i;                                                     
          else if(ioptn==1)                                     
             iadj=i+72;                                                  
             if(iadj>144)
                 iadj=iadj-144;
             end                                                        
          lon=((i-0.5)*2.5)/dlontb((lat)) + 1;
          if(lon>ncells((lat)))
          lon=ncells((lat));     
          end
          eqangl(iadj,j)=arybox(boxnmb(lon,j));
              end
          end
   continue;            %111 CONTINUE
      end
      end
       return;        
                                                                   
%}