#!/bin/ksh 
#
# Abstract: This utility generates a RTG-like SST analysis file with 12Z NSST analysis in FV3GFS.
#
# Notes:
#       RTG-like SST data file: SST analysis over water grids, derived surface temperature over non-water grids.
#                               1/12 lat/lon grids, 4320 x 2160
# Input files :
#       sfcanl      : 6-hourly GFS NSST foundation temperature analysis (tref) and others like surface mask ( land), T1534 Gaussian grids (3072 x 1536)
#       iceanl      : 6-hourly global sea ice analysis, grib (1/12 degree, 4320 x 2160) 
#       grbmsk      : Surface land/water mask 
#       sstclm      : SST climatology (monthly, RTG now)
#       salclm      : Salinity climatology (Monthly, Levitus), at T1534 Gaussian grids. Interpolated from 1 degree lat/lon grids
#                     And then set as 12 for Great Lake, 270 for Sal Lake and 500 for Dead Sea
                                   
# Output file :
#       tf_grb1     : RTG-like Tf analysis with non-water grids filled as RTG, GRIB1
#       tf_grb2     : RTG-like Tf analysis with non-water grids filled as RTG, GRIB2
#
# Author, 03/21/2019     - Xu Li
#


set -x

cd $DATA

msg="HAS BEGUN on `hostname`"
postmsg "$jlogfile" "$msg"

mailbody=mailbody.txt
rm -f $mailbody

#
# SFCIO: ncio or nemsio
#
  SFCIO='ncio'
#
# logical variable if fill land grids or not
#
  lputsi=.true.
#
# search radius
#
  dsearch=500.0
#
# monthly SST climatology (1/12 degree RTG is used here)
#
  SSTCLM=$FIXnsst/RTGSST.1982.2012.monthly.clim.grb
#
# Salinity climatology, half lat/lon degree
#
  SALCLM=$FIXnsst/global_salclm.t1534.3072.1536.nc

  ln -s $SSTCLM    sstclm
  ln -s $SALCLM    salclm
#
# input FV3GFS 6-hourlay surface analysis file at Gaussian grids, for T1534, 3072 x 1536
# If unavailable for the current day, check previous day
#
  for SDIR in $COMINgdas/${cyc} ${COMINgdas/$PDY/$PDYm1}/${cyc}
  do
     if [ $SFCIO = "nemsio" ] ; then
        sfcanl_name=$CDUMP.$cycle.sfcanl.nemsio
     elif [ $SFCIO="ncio" ] ; then
        sfcanl_name=atmos/$CDUMP.$cycle.sfcanl.nc
     else
       echo "invalid sfcio"
       exit
     fi
     if [ -s $SDIR/$sfcanl_name ] ; then
        ln -fs $SDIR/$sfcanl_name sfcanl
        if [ $SDIR != $COMINgdas/${cyc} ]; then 
           echo "**************************************************************************"  > $mailbody
           echo "*** WARNING:  Cannot find $COMINgdas/${cyc}/$sfcanl_name"                   >> $mailbody
           echo "***           Using $SDIR/$sfcanl_name instead"                             >> $mailbody
           echo "**************************************************************************" >> $mailbody
        fi
        break
     fi
  done
  if [ ! -s sfcanl ] ; then
     err_exit "FATAL ERROR: Could not find input file $COMINgdas/${cyc}/$sfcanl_name nor backup $SDIR/$sfcanl_name"
  fi
#
# input Sea ice daily analysis, 1/12 lat/lon grids, 4320 x 2160
# If unavailable for the current day, check previous day
#
  for SDIR in $COMINgdas/${cyc} ${COMINgdas/$PDY/$PDYm1}/${cyc}
  do
     if [ $SFCIO = "nemsio" ] ; then
        iceanl_name=$CDUMP.$cycle.seaice.5min.blend.grb
     elif [ $SFCIO="ncio" ] ; then
        iceanl_name=atmos/$CDUMP.$cycle.seaice.5min.blend.grb
     else
       echo "invalid sfcio for iceanl"
       exit
     fi
     if [ -s $SDIR/$iceanl_name ] ; then
        ln -fs $SDIR/$iceanl_name iceanl
        if [ $SDIR != $COMINgdas/${cyc} ]; then
           echo "**************************************************************************"  > $mailbody
           echo "*** WARNING:  Cannot find $COMINgdas/${cyc}/$iceanl_name"                   >> $mailbody
           echo "***           Using $SDIR/$iceanl_name instead"                             >> $mailbody
           echo "**************************************************************************" >> $mailbody
        fi
        break
     fi
  done
  if [ ! -s iceanl ] ; then
     err_exit "FATAL ERROR: Could not find input file $COMINgdas/${cyc}/$iceanl_name nor backup $SDIR/$iceanl_name"
  fi

# send email if backup data was used
#
  if [ -s "$mailbody" ] ; then
     subject="Missing GDAS data for $PDY t${cyc}z $job"
     mail.py -s "$subject" -c ${MAILCC:-"nco.spa@noaa.gov,Xu.Li@noaa.gov"} < $mailbody
  fi

#
# processing 1/12, 1/4 degree RTG-like or OISST-like file
#
#ResList="twelfth quarter"
 ResList="twelfth quarter"
 for Res in $ResList; do 

    grbmsk=$FIXnsst/rtgssthr_ls_nas.${Res}deg.dat
    ln -fs $grbmsk grbmsk

#
# half degree files are generated while processing 1/12 degree files
#
    if [ $Res = "twelfth" ] ; then
       nx=4320
       ny=2160
       grbmsk_0p5=$FIXnsst/rtgssthr_ls_nas.halfdeg.dat
       ln -fs $grbmsk_0p5 grbmsk_0p5
       ln -fs rtgsstgrb0.083 tf_grb
       ln -fs rtgsstgrb0.083_awips tf_grb_awips
       ln -fs rtgsstgrb0.5 tf_grb_0p5
       ln -fs sstoi_grb tf_grb_1p0
       ln -fs rtgsstgrb0.5_awips tf_grb_0p5_awips
    elif [ $Res = "quarter" ] ; then
       nx=1440
       ny=720
#      ln -fs sstoiqd.sst.grb tf_grb         # old name
       ln -fs nsst.t12z.sst.0p25.grb tf_grb  # new name
    else
       echo "invalid resolution"
       exit
    fi
#
# Make namelist for Res dependent CASE
#
cat << EOF > nsst_parm.input
    &setup
    catime='$CDATE',sfcio='$SFCIO',lputsi=$lputsi,dsearch=$dsearch,nx=$nx,ny=$ny
/
EOF

#
# run nsst executable
#
 $EXECnsst/nsst.x < nsst_parm.input 1>res_${Res}.log 2>res_${Res}.err
 export err=$?;err_chk

 done
 
#
# Get grib1 index files
#
$GRBINDEX rtgsstgrb0.083        rtgsstgrb0.083.index       ;export err=$?;err_chk
$GRBINDEX rtgsstgrb0.083_awips  rtgsstgrb0.083_awips.index ;export err=$?;err_chk
$GRBINDEX rtgsstgrb0.5          rtgsstgrb0.5.index         ;export err=$?;err_chk
$GRBINDEX rtgsstgrb0.5_awips    rtgsstgrb0.5_awips.index   ;export err=$?;err_chk

$GRBINDEX nsst.t12z.sst.0p25.grb nsst.t12z.sst.0p25.grb.index  ;export err=$?;err_chk
$GRBINDEX sstoi_grb              sstoi_grb.index               ;export err=$?;err_chk
#####################################################################
#
# Get grib2 data and index files for awips at 1/12 degree resolution 
#
#####################################################################

$CNVGRIB -g12 -p40 rtgsstgrb0.083_awips  rtgsst_grb_0.083_awips.grib2        ;export err=$?;err_chk
$WGRIB2 rtgsst_grb_0.083_awips.grib2 -s > rtgsst_grb_0.083_awips.grib2.idx   ;export err=$?;err_chk

$CNVGRIB -g12 -p2 rtgsstgrb0.083_awips  rtgsst_grb_0.083_nomads.grib2        ;export err=$?;err_chk
$WGRIB2 rtgsst_grb_0.083_nomads.grib2 -s > rtgsst_grb_0.083_nomads.grib2.idx ;export err=$?;err_chk

export parmcard=${PARMnsst}/grib2_awips_rtgssthr

export pgm=tocgrib
. prep_step
export FORT11="rtgsst_grb_0.083_awips.grib2"
export FORT31=" "
export FORT51="grib2.awips_rtgssthr_grb_0.083"

$TOCGRIB2 < $parmcard
export err=$?;err_chk

#####################################################################
#
# Add WMO Bulletin Header to 1 deg SST field for dissemination
#
#####################################################################
export pgm=tocgrib
. prep_step
export FORT11="sstoi_grb"
export FORT31="sstoi_grb.index"
export FORT51="xtrn.nsst1d_ecmwf"

$TOCGRIB < $PARMnsst/grib_sstecmwf parm='KWBC'
export err=$?;err_chk

######################################################
#
# Convert to grib2 format for 0.5 degree RTG-like file
#
######################################################
$CNVGRIB -g12 -p40 rtgsstgrb0.5 rtgssthr_grb_0.5.grib2              ;export err=$?;err_chk
$CNVGRIB -g12 -p40 rtgsstgrb0.083 rtgssthr_grb_0.083.grib2          ;export err=$?;err_chk
$WGRIB2 rtgssthr_grb_0.5.grib2 -s > rtgssthr_grb_0.5.grib2.idx      ;export err=$?;err_chk
$WGRIB2 rtgssthr_grb_0.083.grib2 -s > rtgssthr_grb_0.083.grib2.idx  ;export err=$?;err_chk
##########################################################
#
# Convert to grib2 format for 0.25 degree OISST-like file
#  - call copygb to create an intermediate grib1 file containing only TMP (11)
#    with land masked out using the landsfc field (81) from the original grib file
#  - convert the masked TMP-only file to grib2 with simple packing (-p0) 
#  - use wgrib2 to change the TMP parameter to WTMP (Water temperature)
#
##########################################################
$COPYGB -x -k'4*-1 11' -B-1 -K'4*-1 81' -A\<0.5 nsst.t12z.sst.0p25.grb nsst.t12z.sst.0p25.grb.masked ;export err=$?;err_chk
$CNVGRIB -g12 -p0 nsst.t12z.sst.0p25.grb.masked nsst.t12z.TMP.0p25.grib2                             ;export err=$?;err_chk
$WGRIB2 nsst.t12z.TMP.0p25.grib2 -if ":TMP:" -set_var WTMP -fi -grib nsst.t12z.sst.0p25.grib2        ;export err=$?;err_chk
$WGRIB2 nsst.t12z.sst.0p25.grib2 -s > nsst.t12z.sst.0p25.grib2.idx                                   ;export err=$?;err_chk
##########################################################
#
# Convert to grib2 format for 1.00 degree OISST-like file
#
##########################################################
$CNVGRIB -g12 -p40 sstoi_grb sstoi_grb.grib2                                           ;export err=$?;err_chk
$WGRIB2 sstoi_grb.grib2 -s > sstoi_grb.grib2.idx                                       ;export err=$?;err_chk
##############################################
#
# Get GEMPAK SST file for atl and pac areas 
#
##############################################
NAGRIB=${GEMEXE}/nagrib_nc
#

  cpyfil=gds
  garea=dset
  gbtbls=
  maxgrd=4999
  kxky=
  grdarea=
  proj=
  output=T

 GRIBIN=rtgsstgrb0.083_awips

 GribList="atl pac"

 for GRIB in $GribList; do

  GEMGRD=rtgsst_${GRIB}_${PDY}00

  if [ $GRIB = "atl" ] ; then
    $COPYGB -g "255 0 840 540 50000 -100000 128 5000 -30000 083 083 0" -x $GRIBIN grib
    export err=$?;err_chk
  else
    # Assume Pacific region
    $COPYGB -g "255 0 960 660 60000 -170000 128 5000 -90000 083 083 0" -x $GRIBIN grib
    export err=$?;err_chk
  fi
  export pgm="nagrib_nc F"
  startmsg

   $NAGRIB << EOF
   GBFILE   = grib
   INDXFL   =
   GDOUTF   = $GEMGRD
   PROJ     = $proj
   GRDAREA  = $grdarea
   KXKY     = $kxky
   MAXGRD   = $maxgrd
   CPYFIL   = $cpyfil
   GAREA    = $garea
   OUTPUT   = $output
   GBTBLS   = $gbtbls
   GBDIAG   =
   PDSEXT   = $pdsext
  l
  r
EOF
  export err=$?;err_chk

done

if test "$SENDCOM" = 'YES'
then
  cp rtgsstgrb0.083_awips             $COMOUT/rtgssthr_grb_0.083_awips
  cp rtgsstgrb0.083_awips.index       $COMOUT/rtgssthr_grb_0.083_awips.index
  cp rtgsstgrb0.083                   $COMOUT/rtgssthr_grb_0.083
  cp rtgsstgrb0.083.index             $COMOUT/rtgssthr_grb_0.083.index
  cp rtgsstgrb0.5_awips               $COMOUT/rtgssthr_grb_0.5_awips
  cp rtgsstgrb0.5_awips.index         $COMOUT/rtgssthr_grb_0.5_awips.index
  cp rtgsstgrb0.5                     $COMOUT/rtgssthr_grb_0.5
  cp rtgsstgrb0.5.index               $COMOUT/rtgssthr_grb_0.5.index

  cp nsst.t12z.sst.0p25.grb           $COMOUT/nsst.t12z.sst.0p25.grb
  cp nsst.t12z.sst.0p25.grb.index     $COMOUT/nsst.t12z.sst.0p25.grb.index
  cp sstoi_grb                        $COMOUT/sstoi_grb
  cp sstoi_grb.index                  $COMOUT/sstoi_grb.index

  cp rtgssthr_grb_0.5.grib2           $COMOUT/rtgssthr_grb_0.5.grib2
  cp rtgssthr_grb_0.083.grib2         $COMOUT/rtgssthr_grb_0.083.grib2
  cp rtgssthr_grb_0.5.grib2.idx       $COMOUT/rtgssthr_grb_0.5.grib2.idx
  cp rtgssthr_grb_0.083.grib2.idx     $COMOUT/rtgssthr_grb_0.083.grib2.idx
  cp rtgsst_grb_0.083_awips.grib2     $COMOUT/rtgssthr_grb_0.083_awips.grib2
  cp rtgsst_grb_0.083_awips.grib2.idx $COMOUT/rtgssthr_grb_0.083_awips.grib2.idx
  cp nsst.t12z.sst.0p25.grib2         $COMOUT/nsst.t12z.sst.0p25.grib2
  cp nsst.t12z.sst.0p25.grib2.idx     $COMOUT/nsst.t12z.sst.0p25.grib2.idx
  cp sstoi_grb.grib2                  $COMOUT/sstoi_grb.grib2
  cp sstoi_grb.grib2.idx              $COMOUT/sstoi_grb.grib2.idx

  cp rtgsst_atl_${PDY}00              $COMOUTgempak/rtgsst_atl_${PDY}00
  cp rtgsst_pac_${PDY}00              $COMOUTgempak/rtgsst_pac_${PDY}00

  cp grib2.awips_rtgssthr_grb_0.083   $COMOUTwmo/grib2.awips_rtgssthr_grb_0.083
  cp xtrn.nsst1d_ecmwf                $COMOUTwmo/xtrn.nsst1d_ecmwf

  if test "$SENDDBN" = 'YES'
  then
     $DBNROOT/bin/dbn_alert MODEL RTG_SST_grib  $job $COMOUT/rtgssthr_grb_0.5
     $DBNROOT/bin/dbn_alert MODEL RTG_SST_gribi $job $COMOUT/rtgssthr_grb_0.5.index
     $DBNROOT/bin/dbn_alert MODEL RTG_SST_grib  $job $COMOUT/rtgssthr_grb_0.083
     $DBNROOT/bin/dbn_alert MODEL RTG_SST_gribi $job $COMOUT/rtgssthr_grb_0.083.index
     $DBNROOT/bin/dbn_alert MODEL RTG_SST_grib  $job $COMOUT/rtgssthr_grb_0.083_awips
     $DBNROOT/bin/dbn_alert MODEL RTG_SST_gribi $job $COMOUT/rtgssthr_grb_0.083_awips.index

     $DBNROOT/bin/dbn_alert MODEL RTGSST_GEMPAK  $job $COMOUTgempak/rtgsst_atl_${PDY}00
     $DBNROOT/bin/dbn_alert MODEL RTGSST_GEMPAK  $job $COMOUTgempak/rtgsst_pac_${PDY}00

#    $DBNROOT/bin/dbn_alert MODEL RTG_SST_grib  $job $COMOUT/sstoi_grb    # reported not needed. commented rather than removed just in case

     $DBNROOT/bin/dbn_alert MODEL RTG_SST_grib_GB2      $job $COMOUT/rtgssthr_grb_0.5.grib2
     $DBNROOT/bin/dbn_alert MODEL RTG_SST_grib_GB2      $job $COMOUT/rtgssthr_grb_0.083.grib2
     $DBNROOT/bin/dbn_alert MODEL RTG_SST_grib_GB2      $job $COMOUT/rtgssthr_grb_0.083_awips.grib2
     $DBNROOT/bin/dbn_alert MODEL RTG_SST_grib_GB2_WIDX $job $COMOUT/rtgssthr_grb_0.5.grib2.idx
     $DBNROOT/bin/dbn_alert MODEL RTG_SST_grib_GB2_WIDX $job $COMOUT/rtgssthr_grb_0.083.grib2.idx
     $DBNROOT/bin/dbn_alert MODEL RTG_SST_grib_GB2_WIDX $job $COMOUT/rtgssthr_grb_0.083_awips.grib2.idx

     $DBNROOT/bin/dbn_alert MODEL NSST_0P25_GB2         $job $COMOUT/nsst.t12z.sst.0p25.grib2

  else
     echo "SENDDBN=$SENDDBN, files not posted to db_net."
  fi

  if test "$SENDDBN_NTC" = 'YES'
  then
     $DBNROOT/bin/dbn_alert NTC_LOW $NET        $job $COMOUTwmo/grib2.awips_rtgssthr_grb_0.083
     $DBNROOT/bin/dbn_alert GRIB_LOW $NET       $job $COMOUTwmo/xtrn.nsst1d_ecmwf
  else
     echo "SENDDBN_NTC=$SENDDBN_NTC, headed files not posted to db_net."
  fi

fi

set -x
msg='ENDED NORMALLY.'
postmsg "$jlogfile" "$msg"

exit
