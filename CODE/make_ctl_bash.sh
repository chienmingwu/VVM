#!/bin/bash

#rundir=`pwd`
rundir=$(cd "$(dirname "$0")" && pwd)   # folder that holds this script
exp=$(echo ${rundir}|rev|cut -d/ -f1|rev)
ncdir="${rundir}/archive"
outdir="${rundir}/gs_ctl_files"
echo ${exp}
echo ${outdir}
mkdir -p ${outdir}

declare -A vtab
vtab[time]='t'
vtab[lon]='x'
vtab[lat]='y'
vtab[lev]='z'

parse_ctl_vars() {
  local fname=$1
  local attr_name=$2
  local outnz=$3
  local type_regex=$4
  local skip_regex=$5

  ncdump -h "${fname}" | awk \
    -v attr_name="${attr_name}" \
    -v outnz="${outnz}" \
    -v type_regex="${type_regex}" \
    -v skip_regex="${skip_regex}" '
    function trim(s) {
      gsub(/^[ \t]+|[ \t]+$/, "", s)
      return s
    }
    function dim_to_ctl(dim, i, n, parts, out, d) {
      gsub(/[ \t]/, "", dim)
      n = split(dim, parts, ",")
      out = ""
      for (i = 1; i <= n; i++) {
        d = parts[i]
        if (d == "time") d = "t"
        else if (d == "lon") d = "x"
        else if (d == "lat") d = "y"
        else if (d == "lev") d = "z"
        out = out (out == "" ? "" : ",") d
      }
      return out
    }
    function quoted_value(line, val) {
      val = line
      sub(/^[^"]*"/, "", val)
      sub(/" ;.*/, "", val)
      return val
    }
    /^[ \t]*(float|int|double)[ \t]+[A-Za-z_][A-Za-z0-9_]*\(/ {
      line = trim($0)
      type = line
      sub(/[ \t].*/, "", type)
      if (type !~ "^(" type_regex ")$") next

      name = line
      sub(/^[^ \t]+[ \t]+/, "", name)
      sub(/\(.*/, "", name)
      if (skip_regex != "" && name ~ "^(" skip_regex ")$") next

      dim = line
      sub(/^[^(]*\(/, "", dim)
      sub(/\).*/, "", dim)

      order[++nvar] = name
      dims[name] = dim_to_ctl(dim)
      longname[name] = name
      next
    }
    index($0, ":" attr_name " = ") {
      line = trim($0)
      name = line
      sub(/:.*/, "", name)
      longname[name] = quoted_value(line)
      next
    }
    /:_FillValue[ \t]*=/ {
      line = trim($0)
      val = line
      sub(/^.*= */, "", val)
      sub(/ ;.*/, "", val)
      sub(/f$/, "", val)
      fillvalue = val
      next
    }
    END {
      if (fillvalue == "") fillvalue = "99999."
      printf "%d %s\n", nvar, fillvalue
      for (i = 1; i <= nvar; i++) {
        name = order[i]
        printf "%s=>%s %s %s %s\n", name, name, outnz, dims[name], longname[name]
      }
    }'
}

# -------------------------------------------
# ----- number of data type and ncheader
dtype_list=""
for dum in $(ls ${ncdir}/*000000.nc);do
  dum=$(echo ${dum}|rev|cut -d"/" -f1|rev|cut -d"-" -f1)
  dum0=$(echo ${dum}|rev|cut -d"." -f2|rev)  #L
  dum1=$(echo ${dum}|rev|cut -d"." -f1|rev)  #Thermodynamic
  ncheader=$(echo ${dum}|rev|cut -d"." -f3|rev)  #RRCE_3km_f00
  dtype_list="${dtype_list} ${dum0}.${dum1}"
done

# -------------------------------------------
# ----- get dimension
nt=$(ls ${ncdir}/*.${dum0}.${dum1}-*.nc|wc -l)
# n=0
# nfmt=$(printf "%06d" ${n})
# echo "${ncdir}/${ncheader}.${dum0}.${dum1}-${nfmt}.nc"
# while [ -f ${ncdir}/${ncheader}.${dum0}.${dum1}-${nfmt}.nc ]; do
#   n=$(echo "${n}+1"|bc)
#   nfmt=$(printf "%06d" ${n})
# done
# nt=${n}

nx=$(grep "zonal_dimension" ${rundir}/DOMAIN|cut -d"'" -f2|cut -d"/" -f3)
ny=$(grep "merid_dimension" ${rundir}/DOMAIN|cut -d"'" -f2|cut -d"/" -f3)
nz=$(grep "vert_dimension" ${rundir}/DOMAIN|cut -d"'" -f2|cut -d"/" -f3)

# -------------------------------------------
# ----- get lon/lat
##  method 1 ---- from INTPUT
### lat0=$(grep "RLAT=" ${rundir}/INPUT|cut -d',' -f1|cut -d'=' -f2)
### lon0=$(grep "RLON=" ${rundir}/INPUT|cut -d',' -f2|cut -d'=' -f2|cut -d' ' -f1)
### #lat0=0
### #lon0=0
### dx=$(grep "DX=" ${rundir}/INPUT|cut -d',' -f1|cut -d'=' -f2)
### dy=$(grep "DY" ${rundir}/INPUT|cut -d',' -f2|cut -d'=' -f2)
### dlon=$(echo "scale=6;${dx}/111000"|bc)
### dlat=$(echo "scale=6;${dy}/111000"|bc)
### 
### lon0=$(echo "scale=3;${lon0}-(${dlon}*${nx}/2)"|bc )
### lat0=$(echo "scale=3;${lat0}-(${dlat}*${ny}/2)"|bc )

### # method 2 ---- from TOPO.nc using ncdump
fname=${ncdir}/${ncheader}.L.Dynamic-000000.nc
lon0=$(ncdump -v,lon ${fname} |grep ' lon ='|cut -c7-|cut -d',' -f1)
lon1=$(ncdump -v,lon ${fname} |grep ' lon ='|cut -c7-|cut -d',' -f2)
lat0=$(ncdump -v,lat ${fname} |grep ' lat ='|cut -c7-|cut -d',' -f1)
lat1=$(ncdump -v,lat ${fname} |grep ' lat ='|cut -c7-|cut -d',' -f2)
dlon=$(echo "scale=7;(${lon1}-(${lon0}))/1"|bc)
dlat=$(echo "scale=7;(${lon1}-(${lon0}))/1"|bc)
lon0=$(printf "%.7f" "$lon0")
lat0=$(printf "%.7f" "$lat0")

# -------------------------------------------
# ----- get dt 
fname=${ncdir}/${ncheader}.L.Dynamic-000001.nc
if [ -f ${fname} ];then
  deltatime=$(ncdump -v,time ${fname} |grep ' time ='|cut -d' ' -f4)
else
  outfreq=$(grep "NXSAVG=" ${rundir}/INPUT |cut -d"," -f2|cut -d"=" -f2)
  dt=$(grep "DT=" ${rundir}/INPUT |cut -d"," -f5|cut -d"=" -f2)
  deltatime=$(echo "scale=0;${outfreq}*${dt}/60"|bc)
fi
deltatime=$(echo "scale=0;${deltatime}/1"|bc)  #keep integer
deltatime=$(( ${deltatime} > 1 ? ${deltatime} : 1 ))  #minmum value is 1
echo "${deltatime} min (${nt})"

# -------------------------------------------
# ----- get z level
# method 1 from fort.98
dum=$(grep -n "ZT(K)" ${rundir}/fort.98|cut -d":" -f1)
dum=$(echo ${dum}+${nz}+1|bc)
table=$(cat ${rundir}/fort.98 |head -n ${dum}|tail -n ${nz})
zlist=""
for i in $(seq ${nz});do
  idx=$(echo "3+(${i}-1)*5"|bc)
  dum=$(echo ${table}|cut -d" " -f${idx})
  dum=$(echo "scale=0;(${dum}/1.)"|bc)
  zlist="${zlist} ${dum}"
  if [ $((i % 10)) -eq 0 ] && [ ${i} -ne ${nz} ]; then
    zlist="${zlist}\n"
  fi
done

#########################################
# create ctl for the data in archive
#########################################
for dtype in ${dtype_list};do
  type0=$(echo ${dtype}|cut -d. -f1)  # L
  type1=$(echo ${dtype}|cut -d. -f2)  # Thermodynamic
  type1_lower=$(echo "${type1}" | sed -e 's/\(.*\)/\L\1/')

  if [ ${type0} == "L" ];then
    outz=${zlist}
    outnz=${nz}
    outnz1=${nz}
  else
    outz=$(echo "${zlist}" | xargs | cut -d' ' -f1)
    outnz=0
    outnz1=1
  fi

  # ------ get varables
  mapfile -t parsed_vars < <(parse_ctl_vars "${ncdir}/${ncheader}.${dtype}-000000.nc" "standard_name" "${outnz}" "float" "time")
  read -r nvar fillvalue <<< "${parsed_vars[0]}"
  varstring=$(printf '\n %s' "${parsed_vars[@]:1}")
  echo ${dtype} ${nvar}

  
  string="
DSET ^../archive/${ncheader}.${dtype}-%tm6.nc
DTYPE netcdf
OPTIONS template
TITLE ${dtype} variables
UNDEF ${fillvalue}
CACHESIZE 10000000
XDEF ${nx} LINEAR ${lon0} ${dlon}
YDEF ${ny} LINEAR ${lat0} ${dlat}
ZDEF ${outnz1} LEVELS ${outz}
TDEF ${nt} LINEAR 01JAN1998 ${deltatime}mn
VARS ${nvar} ${varstring}
ENDVARS
  "
  echo -e "${string}" > "${outdir}/${type1_lower}.ctl"


done


#########################################
# create ctl for TOPO.nc
#########################################
outz=1000
outnz=0
outnz1=1

# ------ get varables
fname="${rundir}/TOPO.nc"
if [ -f ${fname} ]; then
    mapfile -t parsed_vars < <(parse_ctl_vars "${fname}" "long_name" "${outnz}" "float|int|double" "time|lon|lat|lev")
    read -r nvar fillvalue <<< "${parsed_vars[0]}"
    varstring=$(printf '\n %s' "${parsed_vars[@]:1}")
    echo "TOPO.nc ${nvar}"

    string="
DSET ^../TOPO.nc
DTYPE netcdf
TITLE TOPO
UNDEF 99999.
CACHESIZE 10000000
XDEF ${nx} LINEAR ${lon0} ${dlon}
YDEF ${ny} LINEAR ${lat0} ${dlat}
ZDEF ${outnz1} LEVELS ${outz}
TDEF ${nt} LINEAR 01JAN1998 ${deltatime}mn
VARS ${nvar} ${varstring}
ENDVARS
    "
    echo -e "${string}" > "${outdir}/topo.ctl"
else
    echo "skip ... TOPO.nc"
fi

#########################################
## bar.ctl
#########################################
outz=${zlist}
outnz=${nz}
outnz1=${nz}

string="
DSET ^../bar.dat
TITLE mean profile
UNDEF 99999.
XDEF 1 LINEAR ${lon0} ${dlon}
YDEF 1 LINEAR ${lat0} ${dlat}
ZDEF ${outnz1} LEVELS ${outz}
TDEF 1 LINEAR 01JAN1998 ${deltatime}mn
VARS 13 
 pbar   ${outnz} 99 pbar  [Pa]   
 pibar  ${outnz} 99 pibar  
 rho    ${outnz} 99 rho   [kg/m3] 
 th     ${outnz} 99 thbar     [K] 
 qv     ${outnz} 99 qvbar     [kg/kg] 
 UG     ${outnz} 99 UG        [m/s] 
 VG     ${outnz} 99 VG        [m/s] 
 Q1LS   ${outnz} 99 Q1LS      [K/s] 
 Q2LS   ${outnz} 99 Q2LS      [g/g/s] 
 WLS    ${outnz} 99 WLS       [m/s] 
 DZT    ${outnz} 99 delta ZT  [m] 
 the    ${outnz} 99 th_e bar  [K] 
 thes   ${outnz} 99 th_es bar [K] 
ENDVARS
"
echo -e "${string}" > "${outdir}/bar.ctl"
