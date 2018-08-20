#!/usr/bin/perl 

use strict ;
use warnings ;

my %ini_dict_ptr  = ( 'gen_dir' => './gs_ctl_files/' );
initialize( \%ini_dict_ptr ) ;
creat_ctl(\%ini_dict_ptr)    ;
creat_gs(\%ini_dict_ptr)     ;



sub creat_gs{
    my $info = $_[0] ;
    my $color_http = 'http://kodama.fubuki.info/wiki/wiki.cgi/GrADS/script/color.gs?sw=download\&file=color.gs\&count_group=grads' ;
    my $out_name  = $info->{'gen_dir'} . 'color.gs' ;
    print "Getting color script files ...\n" ;
#    `wget -q -O $out_name $color_http` ;

    my %var_nc_hash ;
#    foreach my $cu_nc ( @{$info->{'nclist'}} ){
#        my @var2d = () ;
#        my @var3d = () ;
#        parse_var( \@var2d , \@var3d , $cu_nc ) ;
#    } 

}

sub creat_ctl{
    my $info = $_[0] ;
     
    foreach my $cu_nc ( @{$info->{'nclist'}} ){
        my @var2d = () ;
        my @var3d = () ;
        my $gs_output = ''  ;

        parse_var( \@var2d , \@var3d , $cu_nc ) ;
        
        $gs_output .= 'DSET ^../'.$cu_nc."\n" ;
        $gs_output .= "DTYPE netcdf\n" ;
        $gs_output .= 'TITLE Auto generate ctl by Mars'."\n" ;
        $gs_output .= 'UNDEF 99999.'."\n" ;
        $gs_output .= 'CACHESIZE 10000000'."\n" ;
        $gs_output .= 'XDEF '. $info->{'nx'} . ' LINEAR ' . $info->{'lon'} . ' ' . 
                          scalar( $info->{'dx'} / ( 6378137 * cos( $info->{'lat'} /180 * 3.14159 ) * 2 * 3.14159 ) * 360) .
                          "\n" ;
        $gs_output .= 'YDEF '. $info->{'ny'} . ' LINEAR ' . $info->{'lat'} . ' ' . 
                          scalar ($info->{'dy'} / ( 6378137 * 2 * 3.14159 ) * 360 ) .
                          "\n" ;
        $gs_output .= 'ZDEF '. scalar( @{$info->{'zlevel'}} ) . ' levels ' . "\n"  ;
            foreach ( @{$info->{'zlevel'}} ){
                $gs_output .= $_ . "\n" ;
            }
        $gs_output .= 'TDEF '. scalar( $info->{'ittadd'} / $info->{'nxsavg'}+1 ) . ' LINEAR 00:00Z01May1997 ' .
            scalar( $info->{'nxsavg'} * $info->{'dt'} / 60 ) . "mn \n" ;
        $gs_output .= 'VARS ' . scalar( @var2d + @var3d ) . "\n" ;
        
        foreach ( @var2d  ){
            $gs_output .= $_ . '=>' . $_ . ' 1 t,y,x   From ' . $cu_nc . " \n" ;
        }
        foreach ( @var3d  ){
            $gs_output .= $_ . '=>' . $_ . ' ' .  scalar( @{$info->{'zlevel'}} ) .' t,z,y,x From ' . $cu_nc . " \n" ;
        }
        $gs_output .= 'ENDVARS' ;
           
        if ( $cu_nc =~ m/\.(\w+)\./ ){
            open my $file_ptr , '>' , $info->{'gen_dir'}.$1.".ctl" ;
            print $file_ptr $gs_output ;
        }       
    }
}

sub parse_var{
    my $var2d_list = $_[0] ;
    my $var3d_list = $_[1] ;
    my $nc_name    = $_[2] ;
    
    my $header = `ncdump -h $nc_name` ;
    open my $str_des , '<' , \$header  ;
    while ( <$str_des> ){
    	if (  m/\s+(float|double)\s+(\w+)\((.*)\)\s+;/ ){
            my $dim = $3 ;
            my $varname = $2 ;
            next if ( $varname =~ /(Time|zc|zb|yc|yb|xc|xb)/ ) ;

            if ( $dim =~ m/bottom_top/ ){
                push @$var3d_list , $varname ;
            }else{
                push @$var2d_list , $varname ;
            }
        }
    }
    close $str_des ;
}


sub initialize{
    my $init_dict_ptr = $_[0] ;
    mkdir $init_dict_ptr->{'gen_dir'} ;

    my $raw_dump = `ncdump -v zc -l 10 *.th3d.nc` ;
    if ( $raw_dump =~ m/\s+south_north\s+=\s+(\d+)\s+;/ ){
        $init_dict_ptr->{'ny'} = $1 + 0.0 ;
    }
    if ( $raw_dump =~ m/\s+west_east\s+=\s+(\d+)\s+;/ ){
        $init_dict_ptr->{'nx'} = $1 + 0.0 ;
    }
    $raw_dump =~ s/(.*)zc =(.*);(.*)/$2/gs ;
    $raw_dump =~ s/\s+([0-9.E+]+),*\s*(\n|$)/$1\n/g ;
    my @tempArr = split "\n" , $raw_dump ;
    $init_dict_ptr->{'zlevel'} = \@tempArr ;
 
    my $project_name ;
    my @filelist = split  "\n" , `ls *.th3d.nc` ;
    if ( $filelist[0] =~ m/^(\w+)L\.th3d/ ){
        $init_dict_ptr->{'project_name'} = $1 ;
        $project_name = $1 ;
    }else{
        die "No suitable th3d.nc Files\n" ;
    }
   
    open  my $file_ptr , '<' , $project_name.'.out' ;
    my $flag = 0 ;
    while( <$file_ptr> ){
        if ( m/\s+RLAT\s+=\s+([0-9.E+]+)\s*,?/ ){
            $init_dict_ptr->{'lat'} = $1 + 0.0 ;
            $flag ++ ;
        }
        if ( m/\s+RLON\s+=\s+([0-9.E+]+)\s*,?/ ){
            $init_dict_ptr->{'lon'} = $1 + 0.0 ;
            $flag ++ ;
        }
        if ( m/\s+DX\s+=\s+([0-9.E+]+)\s*,?/ ){
            $init_dict_ptr->{'dx'} = $1 + 0.0 ;
            $flag ++ ;
        }
        if ( m/\s+DYNEW\s+=\s+([0-9.E+]+)\s*,?/ ){
            $init_dict_ptr->{'dy'} = $1 + 0.0 ;
            $flag ++ ;
        }
        if ( m/\s+DT\s+=\s+([0-9.E+]+)\s*,?/ ){
            $init_dict_ptr->{'dt'} = $1 + 0.0 ;
            $flag ++ ;
        }
        if ( m/\s+NXSAVG\s+=\s+([0-9.E+]+)\s*,?/ ){
            $init_dict_ptr->{'nxsavg'} = $1 + 0.0 ;
            $flag ++ ;
        }
        if ( m/\s+ITTADD\s+=\s+([0-9.E+]+)\s*,?/ ){
            $init_dict_ptr->{'ittadd'} = $1 + 0.0 ;
            $flag ++ ;
        }

        last if ( $flag == 7 ) ;
    }
    close $file_ptr ;

    my @nclist = () ;
    push @nclist , $_  while ( <*.nc> ) ;
    $init_dict_ptr->{'nclist'} = \@nclist ;    
    
}
__DATA__

DSET ^../g01C.phys.nc
DTYPE netcdf
TITLE Combine of VVM reslut
UNDEF 99999.
CACHESIZE 10000000
XDEF 128 LINEAR 118.150000 0.018018
YDEF 128 LINEAR 21.250000 0.018018
ZDEF 1   levels 1.000000
TDEF 72 LINEAR 00:00Z01May1997 20mn
VARS 1
sprec=>sprec 1 t,y,x derive from qc3d
ENDVARS

go qc xy contour t
go qv zt grfill t1 t2
