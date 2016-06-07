#!/bin/csh -f 


# global var
setenv master_out "master.dat"	# decide master output name
setenv linewidth  "20"		# set max linewidth
@ icolor=2		# init color index, H.C. recommend to start from red
@ iset=0		# init set index
# seperately legend into box_out
setenv leg_out "legend.dat"   ; set leg_set=0


# check profiles
#set prof_file=(`ls *.prof.*`) ; /bin/echo "found $#prof_file of profiles.."
# set up file name by reaing arguments
if($#argv < 1) then
 cat << eof
 	Usage:
		$0  my_blast_job.out
		
eof
exit
endif

#set prof_file=(`ls $1.prof.*`) ; /bin/echo "found $#prof_file of profiles.."
setenv seed `basename $1 .out`
set prof_file=(`ls $seed.prof.*`) ; /bin/echo "found $#prof_file of profiles.."
set master_prof=(`ls *.prof.0`)
 if($#master_prof == 0) then
   /bin/echo "cant find master profile, continue?(y/n)" 
     setenv opt $< ;  if($opt == "n") exit
   else
    /bin/echo "found master profile: $prof_file[1]"
    set work_prof=`expr $#prof_file - 1` ; /bin/echo "found $work_prof worker profiles"   

   # reset prof_file matrix
   # date: Fri Oct 29 18:50:55 CST 2004
   setenv seed `basename $prof_file[1] .prof.0`
   @ i=1 ; @ j=2
    while($j <= $#prof_file)
      set prof_file[$j]=$seed.prof.$i
     @ i+=1 ; @ j+=1
    end

     @ i=2
      while($i <= $#prof_file)
       @ j=`expr $i - 1`
        /bin/echo "worker # $j : profile -> $prof_file[$i]"
       @ i+=1
      end
 endif
 /bin/echo ""

# dealing with mater and worker seperately
# ------------------------------------------------------------------------
# ------------------		master prof		------------------
# ------------------------------------------------------------------------
# cmd: '{print $2,"1","\n",$2+$(NF),"1","\n&"}'
if(-f $master_out) rm -f $master_out ; setenv jj 0
if(-f $leg_out) rm -f $leg_out

# add pseudo set, to have master stting print out
cat >> $leg_out << eof
# resize legend box, let it fit graph's size
@    view xmin 0.477124
@    view xmax 0.865359
@    view ymin 0.250980
@    view ymax 0.813072
# deactivate the axis, both x and y
@    xaxis  off
@    yaxis  off
# also reset color as white
@    frame color 0
@    s$leg_set line linewidth $linewidth.0
@    s$leg_set line color 0
@    s$leg_set legend  "Master:"
0 $iset
&
eof
@ leg_set+=1

 set rec_result=(`grep "Receive results" $master_prof | awk '{print $2}'`)
 set merg_result=(`grep "Merge fragments" $master_prof | awk '{print $2+$(NF)}'`)
  /bin/echo "found $#rec_result receive-merge patter"
  /bin/echo "dump master profile:" ; /bin/echo ""

# set 0
  cat >> $leg_out << eof
@    s$leg_set line linewidth $linewidth.0
@    s$leg_set line color $icolor
@    s$leg_set legend  "Broadcast Query"
0 $iset
eof
@ leg_set+=1


cat >> $master_out << eof
@    s$iset line linewidth $linewidth.0
@    s$iset line color $icolor
eof
  grep "Broadcast query" $master_prof | awk '{print $2,$2+$(NF)}' | awk '{print $1,jj"\n"$2,jj"\n&"}' jj=$jj >> $master_out

  @ icolor+=1
  @ iset+=1

# multi- set 1
    @ i=1
     while($i <= $#rec_result)
      @ iset+=1
      if($i == 1) then
       cat >> $leg_out << eof 
@    s$leg_set line linewidth $linewidth.0
@    s$leg_set line color $icolor
@    s$leg_set legend  "Rec-Merg Query"
0 $iset
eof
@ leg_set+=1

      endif
       #/bin/echo -n "$rec_result[$i] " ; /bin/echo "scale=20 ; $merg_result[$i] - $rec_result[$i]" | bc
    cat >> $master_out << eof
@    s$iset line linewidth $linewidth.0
@    s$iset line color $icolor
$rec_result[$i] $jj
$merg_result[$i] $jj
&
eof
      @ i+=1
     end  
     @ icolor+=1 
    
    # set 2
    @ iset+=1 
    cat >> $leg_out << eof
@    s$leg_set line linewidth $linewidth.0
@    s$leg_set line color $icolor
@    s$leg_set legend  "Output Results"
0 $iset
eof
@ leg_set+=1

    cat >> $master_out << eof
@    s$iset line linewidth $linewidth.0
@    s$iset line color $icolor
eof

grep "Output results" $master_prof | awk '{print $2,$2+$(NF)}' | awk '{print $1,jj"\n"$2,jj"\n&"}' jj=$jj >> $master_out
 # set up legend ref x position
 set legend_x=`grep "Output results" $master_prof | awk '{print $2}'`

    @ icolor+=1

    # set 3
    @ iset+=1
    cat >> $leg_out << eof
@    s$leg_set line linewidth $linewidth.0
@    s$leg_set line color $icolor
@    s$leg_set legend  "Fetch Outputs"
0 $iset
eof
@ leg_set+=1

    cat >> $master_out << eof
@    s$iset line linewidth $linewidth.0
@    s$iset line color $icolor
eof

	grep "Fetch output"   $master_prof | awk '{print $2,$2+$(NF)}' | awk '{print $1,jj"\n"$2,jj"\n&"}' jj=$jj >> $master_out
     @ icolor+=1

  # fancy plot:
  # feed grace macro definition
#   cp $master_out .tmp
   set nset=(`grep "&" $master_out | wc -l | awk '{print $1}'`)   # grace start from set 0
   echo "total data sets of master profile: $nset" ; @ nset=`expr $nset - 1`
#   rm -f $master_out
#   @ i=1
#    while($i <= $nset)
#     @ j=`expr $i - 1`
#       /bin/echo "@    s$j line linewidth $linewidth.0" >> $master_out
#     @ i+=1
#    end
#    cat .tmp >> $master_out

# add grey pseudo-box to identify worker timeline:

   set box_x=(`grep -B 2 "&" $master_out | grep -v "&" | sed -e 's/--//g' | awk '{print $1}'`)
set y_max=`echo $work_prof + 0.5 | bc`
set y_min="0"
mv $master_out .tmp
 @ i=3 ; @ k=`expr $#box_x - 4`
   while($i < $k)
     @ j=`expr $i + 1`
      set pos_x1=$box_x[$i] ; set pos_x2=$box_x[$j]
     # print only the green portion
cat >> $master_out << eof
@with box
@    box on
@    box loctype world
@    box $pos_x1, $y_max, $pos_x2, $y_min
@    box linestyle 2
@    box linewidth 1.0
@    box color 1
@    box fill color 7
@    box fill pattern 0
@box def
eof
      @ i+=2
   end
   cat .tmp >> $master_out


# ------------------------------------------------------------------------
# ------------------		worker prof		------------------
# ------------------------------------------------------------------------

/bin/echo "" ; /bin/echo "dump worker profile:"
# initial vars:
# reset data set index from max data set of master profile
@ iset=$nset
@ ii=1


if(-f .xmgr) rm -f .xmgr
while($ii <= $work_prof)
@ icolor=2  # reset color index whenever parsing another worker data file
  setenv jj $ii
  # input/output file
  setenv out "worker.$jj.dat" 
  setenv work_prof_fil $seed.prof.$ii
   /bin/echo "" ;/bin/echo "worker # $jj :"; /bin/echo "output file: $out"
   if(-f $out) rm -f $out

   # generate pseudo set for worker, to allow "worker" legend
   if($ii == 1) then
   cat >> $leg_out << eof
@    s$leg_set line linewidth $linewidth.0
@    s$leg_set line color 0
@    s$leg_set legend  "Worker:"
0 $iset
&
eof
@ leg_set+=1

   endif

	# rescale y_pos, for worker prof
	# setenv y_pos `/bin/echo "scale=20; $jj * 0.1" | bc`	
	setenv y_pos $jj

# set 1:
echo "debug: $icolor"
echo "input file: $work_prof_fil"
 set  work_setup=(`grep "File setup took" $work_prof_fil | awk '{print $2}'`)
 # set work_setup_final=(`grep "Fragment list to master took" $work_prof_fil | awk '{print $2+$(NF)}'`) 
 # For temporary use ... because the PROFILE bug in program
 set work_setup_final=(`grep "Fragment list to master took" $work_prof_fil | awk '{print $2}'`) 

# /bin/echo -n "worker setup:" ; /bin/echo "scale=20; $work_setup_final - $work_setup" | bc
#/bin/echo "$work_setup $work_setup_final" > $out
  @ iset+=1
  cat >> $out << eof
@    s$iset line linewidth $linewidth.0
@    s$iset line color $icolor
eof

  if($ii == 1) then 
    cat >> $leg_out << eof
@    s$leg_set line linewidth $linewidth.0
@    s$leg_set line color $icolor
@    s$leg_set legend  "init NCBI "
0 $iset
eof
@ leg_set+=1
  endif

 cat >> $out << eof
$work_setup $y_pos
$work_setup_final $y_pos
&
eof
@ icolor+=1

# set 2
  @ iset+=1 
  cat >> $out << eof
@    s$iset line linewidth $linewidth.0
@    s$iset line color $icolor
eof

# seperate box into another graph
  if($ii == 1) then
    cat >> $leg_out << eof
@    s$leg_set line linewidth $linewidth.0
@    s$leg_set line color $icolor
@    s$leg_set legend  "Rec. Query"
0 $iset
eof
@ leg_set+=1
  endif

  @ icolor+=1

   grep "Receiving query took" $work_prof_fil | awk '{print $2,$2+$(NF)}' | awk '{print $1,jj"\n"$2,jj"\n&"}' jj=$y_pos >> $out

# set 3 & 4
set search_prof=(`grep "Search fragments" $work_prof_fil | awk '{print $2}'`) 
#set cleanup_prof=(`grep "cleanupBLAST() took" $work_prof_fil | awk '{print $2+$(NF)}'`) 
set cleanup_prof=(`grep "cleanupBLAST" $work_prof_fil | awk '{print $2+$(NF)}'`) 

   # print set 3
   set copy_frag=(`grep "Copy fragments" $work_prof_fil | awk '{print $2}'`) 
   set copy_frag_final=(`grep "Copy fragments" $work_prof_fil | awk '{print $2+$(NF)}'`)

 @ i=1
  while($i <= $#copy_frag)
   @ iset+=1

   if($ii == 1 && $i == 1) then
    cat >> $leg_out << eof
@    s$leg_set line linewidth $linewidth.0
@    s$leg_set line color $icolor
@    s$leg_set legend  "Copy Fragment"
0 $iset
eof
@ leg_set+=1
   endif

   cat >> $out << eof
@    s$iset line linewidth $linewidth.0
@    s$iset line color $icolor
$copy_frag[$i] $y_pos
$copy_frag_final[$i] $y_pos
&
eof
    @ i+=1
  end
   @ icolor+=1

  # restart with search-cleanup action, again, with same color set
  @ i=1
  while($i <= $#cleanup_prof)
   @ iset+=1
   cat >> $out << eof
@    s$iset line linewidth $linewidth.0
@    s$iset line color $icolor
eof

    if($ii == 1 && $i == 1) then
     cat >> $leg_out << eof
@    s$leg_set line linewidth $linewidth.0
@    s$leg_set line color $icolor
@    s$leg_set legend  "Search Fragment"
0 $iset
eof
@ leg_set+=1
    endif

   cat >> $out << eof
$search_prof[$i] $y_pos
$cleanup_prof[$i] $y_pos
&
eof
   @ i+=1
  end
   @ icolor+=1

# set 5
# Mon Nov  1 22:17:12 CST 2004 note:
#
# 	Hurng-Chun have update mpiblast-g2 that return to master have been
# 	identified seperately into seach search events.
# jason
# 	update as cascaded portion 

#set final_search=(`grep "Search fragments" $work_prof_fil | tail -1 | awk '{print $2}'`)
#set return_result=(`grep "Return results to master took" $work_prof_fil | awk '{print $2}'`)
#
#  @ iset=`expr $iset + 1` 
#  cat >> $out << eof
#@    s$iset line linewidth $linewidth.0
#@    s$iset line color $icolor
#$final_search $y_pos
#$return_result $y_pos
#&
#eof
# @ icolor+=1
set return_result_ini=(`grep "Return results to master" $work_prof_fil | awk '{print $2}'`) 
set return_result_end=(`grep "Return results to master" $work_prof_fil | awk '{print $2+$(NF)}'`) 
 @ i=1
  while($i <= $#return_result_ini)
   @ iset+=1
 cat >> $out << eof
@    s$iset line linewidth $linewidth.0
@    s$iset line color $icolor
eof

  if($ii == 1 && $i == 1) then
   cat >> $leg_out << eof
@    s$leg_set line linewidth $linewidth.0
@    s$leg_set line color $icolor
@    s$leg_set legend  "Return Results"
0 $iset
eof
@ leg_set+=1
  endif

  cat >> $out << eof
# return result
$return_result_ini[$i] $y_pos
$return_result_end[$i] $y_pos
&
eof
   @ i+=1
  end
  @ icolor+=1


 @ ii+=1

# Thu Oct 28 22:09:56 CST 2004
# skip the following portion, merge into data set, instead of header info
 # fancy plot for worker
 # feed grace macro definition
#  cp $out .tmp
  # accumulate the data set index no
#  @ i=`expr $nset + 1`
  @ itmp=(`grep "&" $out | wc -l | awk '{print $1}'`)
   echo "total no. of data sets: $itmp"
  @ nset=`expr $itmp + $nset`
#   rm -f $out
#    while($i <= $nset)
#     @ j=`expr $i - 1`
#       /bin/echo -e "@    s$j line linewidth $linewidth.0" >> $out
#     @ i+=1
#    end
#    cat .tmp >> $out

 # reset data set index and color index
 @ iset=$nset 

end


# prepare the parameter file
# set x max & min, define as total execution time of master plus tiny step
# and always set xmin as 0.0
set xmax=(`grep -v "&" $master_out | tail -1 | awk '{print $1}'`) ; set xmax=`echo $xmax + 40| bc`
set iwork=(`ls worker.*.dat`) ; set world_max=$#iwork
#set world_max=`echo $#iwork + 0.5 | bc`   # not necessary to def this
#set world_max=`echo $world_max + 4 | bc`  # to fit legend box
 if($world_max >= 16) then
   set world_max=`echo $world_max + 4 | bc`
  else
   set iadd=`echo 16 - $world_max| bc`
   set world_max=`echo $world_max + $iadd + 4 | bc`
 endif

cat >grace.par << eof
world xmin 0.0
world xmax $xmax
world ymin -0.5
world ymax $world_max
xaxis  tick major 600
xaxis  tick minor ticks 9
eof

# ------------------
# | merge all plot |
# ------------------
if(-f plot.dat) rm -f plot.dat
@ spec=`echo $world_max/0.5+2 | bc`
echo "" ;echo "merging master data file..."

#set legend_x=`echo $xmax - 100 | bc # already defiine when parsing master prof
set legend_x=`echo $legend_x + 10 | bc`
# replace determination of legend y as $ymin (=0) + 15, since legend box
# is about 12 in world scale
set legend_y=`echo $world_max - 4 | bc`

mv $leg_out $leg_out.tmp
cat >> $leg_out << eof
# legend box
@    legend on
@    legend loctype view
@    legend 0.5, 0.8
@    legend box color 1
@    legend box pattern 1
@    legend box linewidth 1.0
@    legend box linestyle 1
@    legend box fill color 0
@    legend box fill pattern 1
@    legend font 0
@    legend char size 1.600000
@    legend color 1
@    legend length 2
@    legend vgap 1
@    legend hgap 0
@    legend invert false
eof
cat $leg_out.tmp >> $leg_out  ; rm -f $leg_out.tmp


cat >> plot.dat << eof
# Master data set:
@    xaxis  label "\+\1 Elaspsed Time (sec.)"
#@    world ymin -0.5		# move this to par file
#@    world ymax $world_max	# move this to par file
@    yaxis  tick off
@    yaxis  tick minor off
@    yaxis  tick major 1
@    yaxis  tick spec type both
@    yaxis  tick spec $spec
@    yaxis  tick major 1, 0
@    yaxis  ticklabel 1, "\+\1Master"
eof



 cat $master_out >> plot.dat
  @ i=1
   while($i <= $#iwork)
     echo "merging worker data file # $i ..."
      @ j=`echo "$i / 0.5 + 1" | bc`
     cat >> plot.dat << eof
# worker $i
@    yaxis  tick major $j, $i
@    yaxis  ticklabel $j, "\+\1Worker $i"
eof
     cat worker.$i.dat >> plot.dat
    @ i+=1
   end
    # remove tmp file
    rm $master_out $iwork

    # check if user can access xmgrace 
    # save file as .agr, .pdf, .ps, .png, and raw data "plot.dat"
     which xmgrace > /dev/null
     if($?) then
        /bin/echo -e "cant find command: xmgrace\n exit now\!" ; exit
       else
         #/bin/echo "converting into ps & pdf format..." 
         #xmgrace -hardcopy -printfile plot.ps plot.dat -param grace.par
	 #ps2pdf plot.ps plot.pdf
	 # echo "output files: plot.ps & plot.pdf"
	 #/bin/echo "saving as agr file format: plot.agr" ; xmgrace plot.dat -param grace.par -saveall plot.agr
	 # save as png as well
	 /bin/echo "saving as png file format: plot.png" ; gracebat -hdevice PNG -hardcopy plot.dat -param grace.par
	 echo "saving legend box in other graph: legend.png"
	 gracebat -hdevice PNG -hardcopy $leg_out ; rm $leg_out
         #/bin/echo "starting xmgrace now..." ; xmgrace plot.dat -param grace.par
     endif



# end of script 
# responsible jason
# date: Thu Oct 28 12:23:40 CST 2004
#
# extra notes to access macros of xmgrace:
# 
# 1. max line color is 15
# 2. max line width is 20
# 
# template of master profile
# ------------------------------------------------------------------------
# 0	1.58413	Broadcast query took 0.000556946        (a)  color 1
# 0	35.6684	Receive results from 2 took 0.097404    (b1) c 2
# 0	35.7346	Merge fragments from 2 took 0.000101089 (b1) c 2
# 0	40.5056	Receive results from 1 took 0.127008    (b2)
# 0	40.6014	Merge fragments from 1 took 0.000124931 (b2)
# 0	40.6741	Output results took 0.137194            (c)
# 0	40.8787	Fetch output took 0.376966 		(d)
# 0	41.2558	Total execution time is: 41.2558	(skip)
# ------------------------------------------------------------------------

# template of worker profile
# ------------------------------------------------------------------------
# 1	0.576411	File setup took 0.00920609		<-- ()
# 1	0.59394	initNCBI took 0.0170731				
# 1	0.594396	Fragment list to master took 0.000370074 <--
# 1	0.594425	Receiving query took 0.831334		()
# 1	1.49171	Copy fragments 0 took 38.4541			()
# 1	40.0742	Search fragments 00 				<-- ()
# 1	40.0742	Alias file construction took 0.00025201
# 1	40.1565	InitBlast() took 0.082218
# 1	40.2529	runBLAST() took 0.096268
# 1	40.253	cleanupBLAST() took 3.89769e-05			<--
# 1	40.4122	Search fragments 00				<-- skip
# 1	40.6217	Return results to master took 0.368633		()
# 1	40.6568	File clean up took 0.0350539
# ------------------------------------------------------------------------

# grace color map:
# ------------------------------------------------------------------------
# @map color 0 to (255, 255, 255), "white"
# @map color 1 to (0, 0, 0), "black"
# @map color 2 to (255, 0, 0), "red"
# @map color 3 to (0, 255, 0), "green"
# @map color 4 to (0, 0, 255), "blue"
# @map color 5 to (255, 255, 0), "yellow"
# @map color 6 to (188, 143, 143), "brown"
# @map color 7 to (220, 220, 220), "grey"
# @map color 8 to (148, 0, 211), "violet"
# @map color 9 to (0, 255, 255), "cyan"
# @map color 10 to (255, 0, 255), "magenta"
# @map color 11 to (255, 165, 0), "orange"
# @map color 12 to (114, 33, 188), "indigo"
# @map color 13 to (103, 7, 72), "maroon"
# @map color 14 to (64, 224, 208), "turquoise"
# @map color 15 to (0, 139, 0), "green4"
