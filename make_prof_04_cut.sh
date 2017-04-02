#!/bin/bash
cncfc_dir=/home/adam/Documents/00.projects/02.python/cncfc
#cut files mask
PART=test_1
MOD_SUFIX=2
PART_DIR="part_$PART"
BASE_KNOTS_DIR="$PART_DIR/00.knots"
MODIFIED_KNOTS_DIR="$PART_DIR/01.modif"
GCODE_DIR="$PART_DIR/02.gcode"

BASE_KNOT_PATTERN="????"
PROF_KNOT_PATTERN=????_prof1.knt

KNOT_DOUBLE_1=$BASE_KNOT_PATTERN"d1"
KNOT_DOUBLE_2=$BASE_KNOT_PATTERN"d2"
KNOT_XY=$BASE_KNOT_PATTERN"v?"
KNOT_UV=$BASE_KNOT_PATTERN"u?"
KNOT_FOLD_TOP=$BASE_KNOT_PATTERN"t2"
KNOT_FOLD_SIDE=$BASE_KNOT_PATTERN"s1"
mod_pref=""
cut_xy="xy_"*"$PART""_??.knt"
cut_uv="uv_"*"$PART""_??.knt"
cut_t="$PART"_?t?.knt
cut_c="$PART"_?c?.knt
cut_b="$PART"_?b?.knt
m_cut_t="$mod_pref""$PART"_?t?.knt
m_cut_c="$mod_pref""$PART"_?c?.knt
m_cut_b="$mod_pref""$PART"_?b?.knt

res="$PART*.ngc"
bak="$PART*.bak"
cut="*$PART*$mod_pref*.knt"
res_dir="gcode_"$PART
xyuv_dir="xyuv_"$PART
bak_dir="~bak_"$PART

#if PART dir does not exist, create itss
$cncfc_dir/dxf2knots.py -i test_3.dxf -l prof_04_cutoff
$cncfc_dir/knots2gcode.py -i prof_04_cutoff1.knt -o prof_04_cutoff -d 418 -cm -sh

for FILE in prof_04_xyuv_*.knt
do
  name="${FILE%%.knt}"
  prof_numb="${name##'prof_04_xyuv_'}"
  prof_r="prof_04_r_$prof_numb.knt"
  # echo found matching:
    if [ -f "$prof_r" ]; then
      echo found matching: "$prof_r"

      $cncfc_dir/knots2gcode.py -i $FILE -ir $prof_r -o "prof_04_$prof_numb" -d 418 -cm -sh

    else
      echo did NOT found matching: "$prof_r" - EXIT
      break
    fi
done

cut_file="cut_test4.ngc"
touch $cut_file

cat > $cut_file << EOL
(lofted profiles cut test 1)
G21 G90

F200
G0 B0

G0 X60 Y-5 U60 V-5
EOL

for FILE in prof_04*.ngc
do
echo "o<""${FILE%%.ngc}""> call" >> $cut_file
echo G0 X45 U45 >> $cut_file
echo G0 Y-5 V-5 >> $cut_file
done

cat >> $cut_file << EOL
M2
EOL

# o<prof_02_cutoff> call
