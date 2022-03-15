#!/bin/csh

# ================================================
#  Produce descriptor file for VVM output
# ================================================

# Set plot false if you want to calculate(mass flux, qv flux etc.).
set nonomatch
set plot = true
set dir = archive
set outdir = gs_ctl_files

if( ! -d ${outdir} ) then
  mkdir ${outdir}
endif

if( -f TOPO.nc ) then
  set topo =  true
else
  set topo =  false
endif

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

if( -f *.L.RAS-000000* ) then
  set ras      =  ` ls *.L.RAS-000000* `
else
  set ras      =  nan
endif

if( -f *.L.Diag-000000* ) then
  set diag     =  ` ls *.L.Diag-000000* `
else
  set diag     =  nan
endif

if( -f p3_diagnostic-000000* ) then
  set p3   =  true
else
  set p3   =  false
endif

if( -f surface_diagnostic.dat ) then
  set surdiag   =  true
else
  set surdiag   =  false
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

echo "get dimension informations"

# get information of 3D coordinate  
set dum1    =  ` ncdump -v lon ${dynamic} `
set dum2    =  ` echo ${dum1} | cut -d"=" -f 5 `
set dum3    =  ` echo ${dum1} | cut -d"=" -f 56 `
set nx      =  ` echo ${dum2} | cut -d" " -f 1 `

set dum1    =  ` ncdump -v lat ${dynamic} `
set dum2    =  ` echo ${dum1} | cut -d"=" -f 4 `
set dum3    =  ` echo ${dum1} | cut -d"=" -f 56 `
set ny      =  ` echo ${dum2} | cut -d" " -f 1 `

set dum1    =  ` ncdump -v zc ${dynamic} `
set dum2    =  ` echo ${dum1} | cut -d"=" -f 3 `
set dum3    =  ` echo ${dum1} | rev | cut -d"=" -f 1 | rev `
set nz      =  ` echo ${dum2} | cut -d" " -f 1 `
set prezc   =  ` echo ${dum3} | cut -d";" -f 1 `

set n = 1
set condition = true
set dum3    =  ` echo ${prezc} | cut -d"," -f 1`
set pzc     =  ${dum3}

while ( ${condition} == "true" )
@ n++
set dum3    =  ` echo ${prezc} | cut -d"," -f ${n}`
if( ${dum3} == "" )then
  set condition = false
else
  set pzc   = ` echo ${pzc} ${dum3} `
endif
end

set zc      = ( `echo ${pzc}` )

cd ../${outdir}

cat ../INPUT > input.txt
echo 'PROGRAM xydef \
IMPLICIT NONE \
  \
REAL, PARAMETER :: pi=4*atan(1.) \
REAL :: dx, dy, dz, dz1, lon, lat, cir \
INTEGER :: i, j ,k \
CHARACTER(100) :: var, var2 \
  \
OPEN(10,FILE="input.txt") \
READ(10,*) \
READ(10,999) var2 \
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
i=index(var2,",") \
READ(var2(5:i-1),*) lat \
var2=var2(i+1:100) \
i=index(var2,"/") \
READ(var2(7:i-2),*) lon \
  \
cir=40075.018*cos(lat/180*pi) \
  \
OPEN(10,FILE="xydef.txt") \
WRITE(10,777) lon, ",", dx/1000/cir*360 \
777 FORMAT(f8.4,a1,f12.8) \
WRITE(10,777) lat, ",", dy/1000/40007.86*360 \
CLOSE(10) \
  \
ENDPROGRAM xydef' > xydef.f95 
f95 xydef.f95
./a.out

set prexx   =  ` head -n 1 xydef.txt `
set preyy   =  ` head -n 2 xydef.txt | tail -n 1 `

set xst     =  ` echo ${prexx} | cut -d"," -f 1 `
set xlen    =  ` echo ${prexx} | cut -d"," -f 2 `

set yst     =  ` echo ${preyy} | cut -d"," -f 1 `
set ylen    =  ` echo ${preyy} | cut -d"," -f 2 `

rm a.out input.txt xydef.f95 xydef.txt
cd ../${dir}

echo "get variable informations"

# get information of variables
# =================thermo=========================
set dum1    =  ` ncdump -h ${thermo} | grep float `
set n = 2
set condition = true
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | rev |cut -d" " -f 1 | rev `
set dum4    =  ` echo ${dum1} | cut -d")" -f ${n} | grep -i lev | cut -d"(" -f 1 | rev | cut -d" " -f 1 | rev `
if( ${dum4} == ${dum3} )then
  set prelev = ${nz}
else
  set prelev = '1 '
endif
set prevar  =  ${dum3}

while ( ${condition} == "true" )
@ n++
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | rev |cut -d" " -f 1 | rev `
set dum4    =  ` echo ${dum1} | cut -d")" -f ${n} | grep -i lev | cut -d"(" -f 1 | rev | cut -d" " -f 1 | rev `
if( ${dum3} == ";" )then
  set condition = false
else
  set prevar = ` echo ${prevar} ${dum3} `
endif
if( ${dum4} == ${dum3} )then
  set prelev = ` echo ${prelev} ${nz} `
else
  set prelev = ` echo ${prelev} 1 `
endif
end

set thermovarname = ( `echo ${prevar}` )
set thermovarlev = (`echo ${prelev}` )
# =================dynamic========================
set dum1    =  ` ncdump -h ${dynamic} | grep float `
set n = 2
set condition = true
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | rev |cut -d" " -f 1 | rev `
set dum4    =  ` echo ${dum1} | cut -d")" -f ${n} | grep -i lev | cut -d"(" -f 1 | rev | cut -d" " -f 1 | rev `
if( ${dum4} == ${dum3} )then
  set prelev = ${nz}
else
  set prelev = '1 '
endif
set prevar  =  ${dum3}

while ( ${condition} == "true" )
@ n++
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | rev |cut -d" " -f 1 | rev `
set dum4    =  ` echo ${dum1} | cut -d")" -f ${n} | grep -i lev | cut -d"(" -f 1 | rev | cut -d" " -f 1 | rev `
if( ${dum3} == ";" )then
  set condition = false
else
  set prevar = ` echo ${prevar} ${dum3} `
endif
if( ${dum4} == ${dum3} )then
  set prelev = ` echo ${prelev} ${nz} `
else
  set prelev = ` echo ${prelev} 1 `
endif
end

set dynamicvarname = ( `echo ${prevar}` )
set dynamicvarlev = (`echo ${prelev}` )
# =================surface========================
set dum1    =  ` ncdump -h ${surface} | grep float `
set n = 2
set condition = true
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | rev |cut -d" " -f 1 | rev `
set dum4    =  ` echo ${dum1} | cut -d")" -f ${n} | grep -i lev | cut -d"(" -f 1 | rev | cut -d" " -f 1 | rev `
if( ${dum4} == ${dum3} )then
  set prelev = ${nz}
else
  set prelev = '1 '
endif
set prevar  =  ${dum3}

while ( ${condition} == "true" )
@ n++
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | rev |cut -d" " -f 1 | rev `
set dum4    =  ` echo ${dum1} | cut -d")" -f ${n} | grep -i lev | cut -d"(" -f 1 | rev | cut -d" " -f 1 | rev `
if( ${dum3} == ";" )then
  set condition = false
else
  set prevar = ` echo ${prevar} ${dum3} `
endif
if( ${dum4} == ${dum3} )then
  set prelev = ` echo ${prelev} ${nz} `
else
  set prelev = ` echo ${prelev} 1 `
endif
end

set surfacevarname = ( `echo ${prevar}` )
set surfacevarlev = (`echo ${prelev}` )
# =================land===========================
if( ${land} == "nan" )then
  set landn = ${surface}
else
  set landn = ${land}
endif

set dum1    =  ` ncdump -h ${landn} | grep float `
set n = 2
set condition = true
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | rev |cut -d" " -f 1 | rev `
set dum4    =  ` echo ${dum1} | cut -d")" -f ${n} | grep -i lev | cut -d"(" -f 1 | rev | cut -d" " -f 1 | rev `
if( ${dum4} == ${dum3} )then
  set prelev = ${nz}
else
  set prelev = '1 '
endif
set prevar  =  ${dum3}

while ( ${condition} == "true" )
@ n++
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | rev |cut -d" " -f 1 | rev `
set dum4    =  ` echo ${dum1} | cut -d")" -f ${n} | grep -i lev | cut -d"(" -f 1 | rev | cut -d" " -f 1 | rev `
if( ${dum3} == ";" )then
  set condition = false
else
  set prevar = ` echo ${prevar} ${dum3} `
endif
if( ${dum4} == ${dum3} )then
  set prelev = ` echo ${prelev} ${nz} `
else
  set prelev = ` echo ${prelev} 1 `
endif
end

set landvarname = ( `echo ${prevar}` )
set landvarlev = (`echo ${prelev}` )
# =================radiation======================
if( ${rad} == "nan" )then
  set radn  = ${surface}
else
  set radn  = ${rad}
endif

set dum1    =  ` ncdump -h ${radn} | grep float `
set n = 2
set condition = true
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | rev |cut -d" " -f 1 | rev `
set dum4    =  ` echo ${dum1} | cut -d")" -f ${n} | grep -i lev | cut -d"(" -f 1 | rev | cut -d" " -f 1 | rev `
if( ${dum4} == ${dum3} )then
  set prelev = ${nz}
else
  set prelev = '1 '
endif
set prevar  =  ${dum3}

while ( ${condition} == "true" )
@ n++
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | rev |cut -d" " -f 1 | rev `
set dum4    =  ` echo ${dum1} | cut -d")" -f ${n} | grep -i lev | cut -d"(" -f 1 | rev | cut -d" " -f 1 | rev `
if( ${dum3} == ";" )then
  set condition = false
else
  set prevar = ` echo ${prevar} ${dum3} `
endif
if( ${dum4} == ${dum3} )then
  set prelev = ` echo ${prelev} ${nz} `
else
  set prelev = ` echo ${prelev} 1 `
endif
end

set radvarname = ( `echo ${prevar}` )
set radvarlev = (`echo ${prelev}` )
# =================tracer=========================
if( ${tracer} == "nan" )then
  set tra  = ${surface}
else
  set tra  = ${tracer}
endif

set dum1    =  ` ncdump -h ${tra} | grep float `
set n = 2
set condition = true
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | rev |cut -d" " -f 1 | rev `
set dum4    =  ` echo ${dum1} | cut -d")" -f ${n} | grep -i lev | cut -d"(" -f 1 | rev | cut -d" " -f 1 | rev `
if( ${dum4} == ${dum3} )then
  set prelev = ${nz}
else
  set prelev = '1 '
endif
set prevar  =  ${dum3}

while ( ${condition} == "true" )
@ n++
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | rev |cut -d" " -f 1 | rev `
set dum4    =  ` echo ${dum1} | cut -d")" -f ${n} | grep -i lev | cut -d"(" -f 1 | rev | cut -d" " -f 1 | rev `
if( ${dum3} == ";" )then
  set condition = false
else
  set prevar = ` echo ${prevar} ${dum3} `
endif
if( ${dum4} == ${dum3} )then
  set prelev = ` echo ${prelev} ${nz} `
else
  set prelev = ` echo ${prelev} 1 `
endif
end

set tracervarname = ( `echo ${prevar}` )
set tracervarlev = (`echo ${prelev}` )
# =================RAS============================
if( ${ras} == "nan" )then
  set rasn  = ${surface}
else
  set rasn  = ${ras}
endif

set dum1    =  ` ncdump -h ${rasn} | grep float `
set n = 2
set condition = true
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | rev |cut -d" " -f 1 | rev `
set dum4    =  ` echo ${dum1} | cut -d")" -f ${n} | grep -i lev | cut -d"(" -f 1 | rev | cut -d" " -f 1 | rev `
if( ${dum4} == ${dum3} )then
  set prelev = ${nz}
else
  set prelev = '1 '
endif
set prevar  =  ${dum3}

while ( ${condition} == "true" )
@ n++
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | rev |cut -d" " -f 1 | rev `
set dum4    =  ` echo ${dum1} | cut -d")" -f ${n} | grep -i lev | cut -d"(" -f 1 | rev | cut -d" " -f 1 | rev `
if( ${dum3} == ";" )then
  set condition = false
else
  set prevar = ` echo ${prevar} ${dum3} `
endif
if( ${dum4} == ${dum3} )then
  set prelev = ` echo ${prelev} ${nz} `
else
  set prelev = ` echo ${prelev} 1 `
endif
end

set rasvarname = ( `echo ${prevar}` )
set rasvarlev = (`echo ${prelev}` )
# =================DIAG===========================
if( ${diag} == "nan" )then
  set diagn  = ${surface}
else
  set diagn  = ${diag}
endif

set dum1    =  ` ncdump -h ${diagn} | grep float `
set n = 2
set condition = true
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | rev |cut -d" " -f 1 | rev `
set dum4    =  ` echo ${dum1} | cut -d")" -f ${n} | grep -i lev | cut -d"(" -f 1 | rev | cut -d" " -f 1 | rev `
if( ${dum4} == ${dum3} )then
  set prelev = ${nz}
else
  set prelev = '1 '
endif
set prevar  =  ${dum3}

while ( ${condition} == "true" )
@ n++
set dum2    =  ` echo ${dum1} | cut -d"(" -f ${n} `
set dum3    =  ` echo ${dum2} | rev |cut -d" " -f 1 | rev `
set dum4    =  ` echo ${dum1} | cut -d")" -f ${n} | grep -i lev | cut -d"(" -f 1 | rev | cut -d" " -f 1 | rev `
if( ${dum3} == ";" )then
  set condition = false
else
  set prevar = ` echo ${prevar} ${dum3} `
endif
if( ${dum4} == ${dum3} )then
  set prelev = ` echo ${prelev} ${nz} `
else
  set prelev = ` echo ${prelev} 1 `
endif
end

set diagvarname = ( `echo ${prevar}` )
set diagvarlev = (`echo ${prelev}` )
# ================================================


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
ZDEF '${nz}' LEVELS ' > thermodynamic.ctl

set n = 0
while ( ${n} < $#zc )
@ n++
echo ${zc[${n}]} >> thermodynamic.ctl
end

echo 'TDEF '${nt}' LINEAR 00:00Z01JAN2000 1mn \
VARS '$#thermovarname >> thermodynamic.ctl

set n = 0
while ( ${n} < $#thermovarname )
@ n++
echo ${thermovarname[${n}]}'=>'${thermovarname[${n}]}' '${thermovarlev[${n}]}' t,z,y,x '${expname} >> thermodynamic.ctl
end
echo 'ENDVARS' >> thermodynamic.ctl
# =================p3_diagnostic==================================
if( ${p3} == "true" )then
echo 'DSET ^../archive/p3_diagnostic-%tm6.dat \
OPTIONS template \
TITLE thermodynamic variables \
UNDEF -99999. \
XDEF '${nx}' LINEAR '${xst}' '${xlen}' \
YDEF '${ny}' LINEAR '${yst}' '${ylen}' \
ZDEF '${nz}' LEVELS ' > p3_diagnostic.ctl

set n = 0
while ( ${n} < $#zc )
@ n++
echo ${zc[${n}]} >> p3_diagnostic.ctl
end

echo 'TDEF '${nt}' LINEAR 00:00Z01JAN2000 1mn \
VARS 5 \
vmi  '${nz}' 99 mass-weighted fall speed (ice) ms-1 \
effi '${nz}' 99 effective radius (ice) m \
di   '${nz}' 99 mean diameter (ice) m \
rhoi '${nz}' 99 bulk density (ice) kgm-1 \
ze   '${nz}' 99 equivalent reflectivity dBz \
ENDVARS' >> p3_diagnostic.ctl
endif
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
ZDEF '${nz}' LEVELS ' > dynamic.ctl

set n = 0
while ( ${n} < $#zc )
@ n++
echo ${zc[${n}]} >> dynamic.ctl
end

echo 'TDEF '${nt}' LINEAR 00:00Z01JAN2000 1mn \
VARS '$#dynamicvarname >> dynamic.ctl

set n = 0
while ( ${n} < $#dynamicvarname )
@ n++
echo ${dynamicvarname[${n}]}'=>'${dynamicvarname[${n}]}' '${dynamicvarlev[${n}]}' t,z,y,x '${expname} >> dynamic.ctl
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
echo ${surfacevarname[${n}]}'=>'${surfacevarname[${n}]}' '${surfacevarlev[${n}]}' t,y,x '${expname} >> surface.ctl
end
echo 'ENDVARS' >> surface.ctl
# =================Surface==================================
if( ${surdiag} == "true" )then
echo 'DSET ^../archive/surface_diagnostic.dat \
TITLE surface variables \
UNDEF -99999 \
XDEF '${nx}' LINEAR '${xst}' '${xlen}' \
YDEF '${ny}' LINEAR '${yst}' '${ylen}' \
ZDEF 1 LEVELS 0 \
TDEF '${nt}' LINEAR 00:00Z01JAN2000 1mn \
VARS 4 \
th  1  99 surface potential temperature (K) \
u   1  99 surface zonal wind (m/s) \
v   1  99 surface meridional wind (m/s) \
w   1  99 surface vertical wind (m/s)\
ENDVARS' > sur_diag.ctl
endif
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
echo ${landvarname[${n}]}'=>'${landvarname[${n}]}' '${landvarlev[${n}]}' t,y,x '${expname} >> land.ctl
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
ZDEF '${nz}' LEVELS ' > radiation.ctl

set n = 0
while ( ${n} < $#zc )
@ n++
echo ${zc[${n}]} >> radiation.ctl
end

echo 'TDEF '${nt}' LINEAR 00:00Z01JAN2000 1mn \
VARS '$#radvarname >> radiation.ctl

set n = 0
while ( ${n} < $#radvarname )
@ n++
echo ${radvarname[${n}]}'=>'${radvarname[${n}]}' '${radvarlev[${n}]}' t,z,y,x '${expname} >> radiation.ctl
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
ZDEF '${nz}' LEVELS ' > tracer.ctl

set n = 0
while ( ${n} < $#zc )
@ n++
echo ${zc[${n}]} >> tracer.ctl
end

echo 'TDEF '${nt}' LINEAR 00:00Z01JAN2000 1mn \
VARS '$#tracervarname >> tracer.ctl

set n = 0
while ( ${n} < $#tracervarname )
@ n++
echo ${tracervarname[${n}]}'=>'${tracervarname[${n}]}' '${tracervarlev[${n}]}' t,z,y,x '${expname} >> tracer.ctl
end
echo 'ENDVARS' >> tracer.ctl
endif
# =================RAS======================================
if( ${ras} != "nan" )then
echo 'DSET ^../archive/'${expname}'.L.RAS-%tm6'${tail}' \
DTYPE netcdf \
OPTIONS template \
TITLE tracers \
UNDEF 9.96921e+36 \
CACHESIZE 10000000 \
XDEF '${nx}' LINEAR '${xst}' '${xlen}' \
YDEF '${ny}' LINEAR '${yst}' '${ylen}' \
ZDEF '${nz}' LEVELS ' > ras.ctl

set n = 0
while ( ${n} < $#zc )
@ n++
echo ${zc[${n}]} >> ras.ctl
end

echo 'TDEF '${nt}' LINEAR 00:00Z01JAN2000 1mn \
VARS '$#rasvarname >> ras.ctl

set n = 0
while ( ${n} < $#rasvarname )
@ n++
echo ${rasvarname[${n}]}'=>'${rasvarname[${n}]}' '${rasvarlev[${n}]}' t,z,y,x '${expname} >> ras.ctl
end
echo 'ENDVARS' >> ras.ctl
endif
# =================RAS======================================
if( ${diag} != "nan" )then
echo 'DSET ^../archive/'${expname}'.L.Diag-%tm6'${tail}' \
DTYPE netcdf \
OPTIONS template \
TITLE tracers \
UNDEF 9.96921e+36 \
CACHESIZE 10000000 \
XDEF '${nx}' LINEAR '${xst}' '${xlen}' \
YDEF '${ny}' LINEAR '${yst}' '${ylen}' \
ZDEF '${nz}' LEVELS ' > diag.ctl

set n = 0
while ( ${n} < $#zc )
@ n++
echo ${zc[${n}]} >> diag.ctl
end

echo 'TDEF '${nt}' LINEAR 00:00Z01JAN2000 1mn \
VARS '$#diagvarname >> diag.ctl

set n = 0
while ( ${n} < $#diagvarname )
@ n++
echo ${diagvarname[${n}]}'=>'${diagvarname[${n}]}' '${diagvarlev[${n}]}' t,z,y,x '${expname} >> diag.ctl
end
echo 'ENDVARS' >> diag.ctl
endif
# =================topography===============================
if( ${topo} == "true" )then
echo 'DSET ^../TOPO.nc \
DTYPE netcdf \
TITLE TOPOGRAPHY \
UNDEF -99999. \
CACHESIZE 10000000 \
XDEF '${nx}' LINEAR '${xst}' '${xlen}' \
YDEF '${ny}' LINEAR '${yst}' '${ylen}' \
ZDEF 1 LEVELS 0 \
TDEF 1 LINEAR 00:00Z01JAN2000 1mn \
VARS 9 \
TOPO=>topo 1 y,x topo \
albedo=>albedo 1 y,x topo \
GRF=>grf 1 y,x topo \
LAI=>lai 1 y,x topo \
LU=>lu 1 y,x topo \
SHDMAX=>shdmax 1 y,x topo \
SHDMIN=>shdmin 1 y,x topo \
SLOPE=>slope 1 y,x topo\
SOIL=>soil 1 y,x topo \
ENDVARS' > topo.ctl
endif
# ==========================================================

exit
