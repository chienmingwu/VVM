#!/bin/csh

# ================================================
#  Produce descriptor file for VVM output
# ================================================

# Set plot false if you want to calculate(mass flux, qv flux etc.).
set plot = true
set dir = archive
set outdir = gs_ctl_files

if( -d ${outdir} ) then
  rm -r ${outdir}
endif
mkdir ${outdir}
cd ${dir}

set   thermo   =  ` ls *.L.Thermodynamic-000000* `
set   dynamic  =  ` ls *.L.Dynamic-000000* `
set   surface  =  ` ls *.C.Surface-000000* `

if( -f *.C.LandSurface-000000* ) then
  set land     =  ` ls *.C.LandSurface-000000* `
else
  set land     =  nan
endif

if( -f *.L.Radiation-000000* ) then
  set rad      =  ` ls *.L.Radiation-000000* `
else
  set rad      =  nan
endif

if( -f *.L.Tracer-000000* ) then
  set tracer   =  ` ls *.L.Tracer-000000* `
else
  set tracer   =  nan
endif

set tail       =  ` echo ${thermo} | rev | cut -d"0" -f 1 | rev `

if( ${plot} == "true" )then
  set gridtype =  agrid 
else
  set gridtype =  ` echo ${thermo} | cut -d"." -f 4 `
endif
 
# define final output file's name
set n = 0
set condition = true

while ( ${condition} == "true" )
@ n++
set number = ` printf "%06d" ${n} `
if( ! -f *.L.Thermodynamic-${number}* )then
  set condition = false
  set nt = ${n}
endif
end

# get information of 3D coordinate  
set dum1    =  ` ncdump -v xc ${dynamic} `
set dum2    =  ` echo ${dum1} | cut -d"=" -f 5 `
set dum3    =  ` echo ${dum1} | cut -d"=" -f 56 `
set nx      =  ` echo ${dum2} | cut -d" " -f 1 `
set dum4    =  ` echo ${dum3} | cut -d";" -f 1 `
set dum5    =  ` echo ${dum4} | cut -d"," -f 1 `
set dum6    =  ` echo ${dum5} | rev | cut -c 1-3 | rev `
set dum7    =  ` echo ${dum5} | rev | cut -c 4-10 | rev `
if ( ${dum7} < 1 )then
set dum7    = 0
endif
@   dum8    = ${dum5} * 2
set dum9    =  ` echo ${dum8} | rev | cut -c 1-3 | rev `
set dum10   =  ` echo ${dum8} | rev | cut -c 4-10 | rev `
if ( ${dum10} < 1 )then
set dum10   = 0
endif

set xst     =  ${dum7}'.'${dum6} 
set xlen    =  ${dum10}'.'${dum9}

set dum1    =  ` ncdump -v yc ${dynamic} `
set dum2    =  ` echo ${dum1} | cut -d"=" -f 4 `
set dum3    =  ` echo ${dum1} | cut -d"=" -f 56 `
set ny      =  ` echo ${dum2} | cut -d" " -f 1 `
set dum4    =  ` echo ${dum3} | cut -d";" -f 1 `
set dum5    =  ` echo ${dum4} | cut -d"," -f 1 `
set dum6    =  ` echo ${dum5} | rev | cut -c 1-3 | rev `
set dum7    =  ` echo ${dum5} | rev | cut -c 4-10 | rev `
if ( ${dum7} < 1 )then
set dum7    = 0
endif
@   dum8    = ${dum5} * 2
set dum9    =  ` echo ${dum8} | rev | cut -c 1-3 | rev `
set dum10   =  ` echo ${dum8} | rev | cut -c 4-10 | rev `
if ( ${dum10} < 1 )then
set dum10   = 0
endif

set yst     =  ${dum7}'.'${dum6}
set ylen    =  ${dum10}'.'${dum9}

set dum1    =  ` ncdump -v zc ${dynamic} `
set dum2    =  ` echo ${dum1} | cut -d"=" -f 3 `
set dum3    =  ` echo ${dum1} | cut -d"=" -f 56 `
set nz      =  ` echo ${dum2} | cut -d" " -f 1 `
set zc      =  ` echo ${dum3} | cut -d";" -f 1 `


# get information of variables
# =================thermo=========================
set dum1    =  ` ncdump -h ${thermo} | grep float `
set n = 2
set condition = true
set dum2    =  ` echo ${dum1} | cut -d"(" -f 2 `
set dum3    =  ` echo ${dum2} | cut -d" " -f 4 `
set prevar  =  ${dum3}

while ( ${condition} == "true" )
@ n++
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | cut -d" " -f 7 `
if( ${dum3} == "" )then
  set condition = false
else
  set prevar = ` echo ${prevar} ${dum3} `
endif
end

set thermovarname = ( `echo ${prevar}` )
# =================dynamic========================
set dum1    =  ` ncdump -h ${dynamic} | grep float `
set n = 2
set condition = true
set dum2    =  ` echo ${dum1} | cut -d"(" -f 2 `
set dum3    =  ` echo ${dum2} | cut -d" " -f 4 `
set prevar  =  ${dum3}

while ( ${condition} == "true" )
@ n++
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | cut -d" " -f 7 `
if( ${dum3} == "" )then
  set condition = false
else
  set prevar = ` echo ${prevar} ${dum3} `
endif
end

set dynamicvarname = ( `echo ${prevar}` )
# =================surface========================
set dum1    =  ` ncdump -h ${surface} | grep float `
set n = 2
set condition = true
set dum2    =  ` echo ${dum1} | cut -d"(" -f 2 `
set dum3    =  ` echo ${dum2} | cut -d" " -f 4 `
set prevar  =  ${dum3}

while ( ${condition} == "true" )
@ n++
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | cut -d" " -f 6 `
if( ${dum3} == "" )then
  set condition = false
else
  set prevar = ` echo ${prevar} ${dum3} `
endif
end

set surfacevarname = ( `echo ${prevar}` )
# =================land===========================
if( ${land} != "nan" )then
set dum1    =  ` ncdump -h ${land} | grep float `
set n = 2
set condition = true
set dum2    =  ` echo ${dum1} | cut -d"(" -f 2 `
set dum3    =  ` echo ${dum2} | cut -d" " -f 4 `
set prevar  =  ${dum3}

while ( ${condition} == "true" )
@ n++
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | cut -d" " -f 6 `
if( ${dum3} == "" )then
  set condition = false
else
  set prevar = ` echo ${prevar} ${dum3} `
endif
end

set landvarname = ( `echo ${prevar}` )
endif
# =================radiation======================
if( ${rad} != "nan" )then
set dum1    =  ` ncdump -h ${rad} | grep float `
set n = 2
set condition = true
set dum2    =  ` echo ${dum1} | cut -d"(" -f 2 `
set dum3    =  ` echo ${dum2} | cut -d" " -f 4 `
set prevar  =  ${dum3}

while ( ${condition} == "true" )
@ n++
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | cut -d" " -f 7 `
if( ${dum3} == "" )then
  set condition = false
else
  set prevar = ` echo ${prevar} ${dum3} `
endif
end

set radvarname = ( `echo ${prevar}` )
endif
# =================tracer=========================
if( ${tracer} != "nan" )then
set dum1    =  ` ncdump -h ${tracer} | grep float `
set n = 2
set condition = true
set dum2    =  ` echo ${dum1} | cut -d"(" -f 2 `
set dum3    =  ` echo ${dum2} | cut -d" " -f 4 `
set prevar  =  ${dum3}

while ( ${condition} == "true" )
@ n++
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | cut -d" " -f 7 `
if( ${dum3} == "" )then
  set condition = false
else
  set prevar = ` echo ${prevar} ${dum3} `
endif
end

set tracervarname = ( `echo ${prevar}` )
endif



cd ../${outdir}
set expname = ` echo ${thermo} | cut -d"." -f 1 `

# produce ctl 
# =================Thermodynamic============================
echo 'DSET ^../archive/'${expname}'.L.Thermodynamic-%tm6'${tail}' \
DTYPE netcdf \
OPTIONS template \
TITLE thermodynamic variables \
UNDEF 99999. \
CACHESIZE 10000000 \
XDEF '${nx}' LINEAR '${xst}' '${xlen}' \
YDEF '${ny}' LINEAR '${yst}' '${ylen}' \
ZDEF '${nz}' LEVELS '${zc}' \
TDEF '${nt}' LINEAR 00:00Z01JAN2000 1mn \
VARS '$#thermovarname > thermodynamic.ctl

set n = 0
while ( ${n} < $#thermovarname )
@ n++
echo ${thermovarname[${n}]}'=>'${thermovarname[${n}]}' '${nz}' t,z,y,x '${expname} >> thermodynamic.ctl
end
echo 'ENDVARS' >> thermodynamic.ctl
# =================Dynamic==================================
if( ${gridtype} == "agrid" )then
# A-grid====================================================
echo 'DSET ^../archive/'${expname}'.L.Dynamic-%tm6'${tail}' \
DTYPE netcdf \
OPTIONS template \
TITLE dynamic variables \
UNDEF 99999. \
CACHESIZE 10000000 \
XDEF '${nx}' LINEAR '${xst}' '${xlen}' \
YDEF '${ny}' LINEAR '${yst}' '${ylen}' \
ZDEF '${nz}' LEVELS '${zc}' \
TDEF '${nt}' LINEAR 00:00Z01JAN2000 1mn \
VARS '$#dynamicvarname > dynamic.ctl

set n = 0
while ( ${n} < $#dynamicvarname )
@ n++
echo ${dynamicvarname[${n}]}'=>'${dynamicvarname[${n}]}' '${nz}' t,z,y,x '${expname} >> dynamic.ctl
end
echo 'ENDVARS' >> dynamic.ctl
else
# C-grid====================================================
cat ../INPUT > input.txt
echo 'PROGRAM cgrid \
IMPLICIT NONE \
  \
REAL*8 :: dx, dy, dz, dz1, cz1, cz2, zb, domain \
INTEGER :: i, j ,k \
REAL*8, DIMENSION('${nz}') :: x, y, zc, ze \
CHARACTER(100) :: var \
  \
OPEN(10,FILE="input.txt") \
READ(10,*) \
READ(10,*) \
READ(10,999) var \
999 FORMAT(10x,100A) \
CLOSE(10) \
  \
i=index(var,",") \
READ(var(3:i-1),*) dx \
var=var(i+1:100) \
i=index(var,",") \
READ(var(7:i-1),*) dy \
var=var(i+1:100) \
i=index(var,",") \
READ(var(4:i-1),*) dz \
var=var(i+1:100) \
i=index(var,",") \
READ(var(5:i-1),*) dz1 \
  \
zb=0. \
domain=15000.  \
cz2=(dz-dz1)/(dz*(domain-dz)) \
cz1=1.-cz2*domain \
  \
ze(1) = 0 \
DO i=2,'${nz}' \
  ze(i) = ze(i-1) + dz \
ENDDO \
  \
zc(1) = ze(1) \
zc(2) = ze(1) + dz*0.5 \
DO i=3,'${nz}' \
  zc(i) = zc(i-1) + dz \
ENDDO \
  \
DO i=1,'${nz}' \
  ze(i) = ze(i) * (cz1+cz2*ze(i)) \
  zc(i) = zc(i) * (cz1+cz2*zc(i)) \
ENDDO \
  \
OPEN(10,FILE="xi.txt") \
WRITE(10,*) dx/1000/2, dx/1000 \
WRITE(10,*) dy/1000, dy/1000 \
WRITE(10,*) ze \
CLOSE(10) \
  \
OPEN(10,FILE="eta.txt") \
WRITE(10,*) dx/1000, dx/1000 \
WRITE(10,*) dy/1000/2, dy/1000 \
WRITE(10,*) ze \
CLOSE(10) \
  \
OPEN(10,FILE="zeta.txt") \
WRITE(10,*) dx/1000, dx/1000 \
WRITE(10,*) dy/1000, dy/1000 \
WRITE(10,*) zc \
CLOSE(10) \
  \
OPEN(10,FILE="u.txt") \
WRITE(10,*) dx/1000, dx/1000 \
WRITE(10,*) dy/1000/2, dy/1000 \
WRITE(10,*) zc \
CLOSE(10) \
  \
OPEN(10,FILE="v.txt") \
WRITE(10,*) dx/1000/2, dx/1000 \
WRITE(10,*) dy/1000, dy/1000 \
WRITE(10,*) zc \
CLOSE(10) \
  \
OPEN(10,FILE="w.txt") \
WRITE(10,*) dx/1000/2, dx/1000 \
WRITE(10,*) dy/1000/2, dy/1000 \
WRITE(10,*) ze \
CLOSE(10) \
  \
ENDPROGRAM cgrid' > cgrid.f95  
f95 cgrid.f95
./a.out

set filename = ("zdx" "zdy" "zdz" "u" "v" "w")
set varname = ("xi" "eta" "zeta" "u" "v" "w")
set n = 1
while ( ${n} <= 6 )
set xxx = ` head -n 1 ${varname[${n}]}'.txt' `
set yyy = ` head -n 2 ${varname[${n}]}'.txt' | tail -n 1 `
set zzz = ` head -n 3 ${varname[${n}]}'.txt' | tail -n 1 `

echo 'DSET ^../archive/'${expname}'.L.Dynamic-%tm6'${tail}' \
DTYPE netcdf \
OPTIONS template \
TITLE dynamic variables \
UNDEF 9.96921e+36 \
CACHESIZE 10000000 \
XDEF '${nx}' LINEAR '${xxx}' \
YDEF '${ny}' LINEAR '${yyy}' \
ZDEF '${nz}' LEVELS '${zzz}' \
TDEF '${nt}' LINEAR 00:00Z01JAN2000 1mn \
VARS  1' > ${filename[${n}]}.ctl
echo ${varname[${n}]}'=>'${varname[${n}]}' '${nz}' t,z,y,x '${expname} >> ${filename[${n}]}.ctl
echo 'ENDVARS' >> ${filename[${n}]}.ctl
@ n++
end
rm *.txt a.out cgrid.f95
endif
# =================Surface==================================
echo 'DSET ^../archive/'${expname}'.C.Surface-%tm6'${tail}' \
DTYPE netcdf \
OPTIONS template \
TITLE surface variables \
UNDEF 9.96921e+36 \
CACHESIZE 10000000 \
XDEF '${nx}' LINEAR '${xst}' '${xlen}' \
YDEF '${ny}' LINEAR '${yst}' '${ylen}' \
ZDEF 1 LEVELS 0 \
TDEF '${nt}' LINEAR 00:00Z01JAN2000 1mn \
VARS '$#surfacevarname > surface.ctl

set n = 0
while ( ${n} < $#surfacevarname )
@ n++
echo ${surfacevarname[${n}]}'=>'${surfacevarname[${n}]}' 1 t,y,x '${expname} >> surface.ctl
end
echo 'ENDVARS' >> surface.ctl
# =================Land=====================================
if( ${land} != "nan" )then
echo 'DSET ^../archive/'${expname}'.C.LandSurface-%tm6'${tail}' \
DTYPE netcdf \
OPTIONS template \
TITLE land surface variables \
UNDEF 9.96921e+36 \
CACHESIZE 10000000 \
XDEF '${nx}' LINEAR '${xst}' '${xlen}' \
YDEF '${ny}' LINEAR '${yst}' '${ylen}' \
ZDEF 1 LEVELS 0 \
TDEF '${nt}' LINEAR 00:00Z01JAN2000 1mn \
VARS '$#landvarname > land.ctl

set n = 0
while ( ${n} < $#landvarname )
@ n++
echo ${landvarname[${n}]}'=>'${landvarname[${n}]}' 1 t,y,x '${expname} >> land.ctl
end
echo 'ENDVARS' >> land.ctl
endif
# =================Radiation================================
if( ${rad} != "nan" )then
echo 'DSET ^../archive/'${expname}'.L.Radiation-%tm6'${tail}' \
DTYPE netcdf \
OPTIONS template \
TITLE radiation variables \
UNDEF 9.96921e+36 \
CACHESIZE 10000000 \
XDEF '${nx}' LINEAR '${xst}' '${xlen}' \
YDEF '${ny}' LINEAR '${yst}' '${ylen}' \
ZDEF '${nz}' LEVELS '${zc}' \
TDEF '${nt}' LINEAR 00:00Z01JAN2000 1mn \
VARS '$#radvarname > radiation.ctl

set n = 0
while ( ${n} < $#radvarname )
@ n++
echo ${radvarname[${n}]}'=>'${radvarname[${n}]}' '${nz}' t,z,y,x '${expname} >> radiation.ctl
end
echo 'ENDVARS' >> radiation.ctl
endif
# =================Tracer===================================
if( ${tracer} != "nan" )then
echo 'DSET ^../archive/'${expname}'.L.Tracer-%tm6'${tail}' \
DTYPE netcdf \
OPTIONS template \
TITLE tracers \
UNDEF 9.96921e+36 \
CACHESIZE 10000000 \
XDEF '${nx}' LINEAR '${xst}' '${xlen}' \
YDEF '${ny}' LINEAR '${yst}' '${ylen}' \
ZDEF '${nz}' LEVELS '${zc}' \
TDEF '${nt}' LINEAR 00:00Z01JAN2000 1mn \
VARS '$#tracervarname > tracer.ctl

set n = 0
while ( ${n} < $#tracervarname )
@ n++
echo ${tracervarname[${n}]}'=>'${tracervarname[${n}]}' '${nz}' t,z,y,x '${expname} >> tracer.ctl
end
echo 'ENDVARS' >> tracer.ctl
endif
# ==========================================================

exit
