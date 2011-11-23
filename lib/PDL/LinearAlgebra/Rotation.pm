package PDL::LinearAlgebra::Rotation;
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;

our $PI=2*acos(0);

sub rotation_matrix_axis_angle {
	my ($u,$angle)=@_;
	return(identity(3)) if (abs($angle) < 0.0000001);
	# shamelessly stolen from wikipedia
	my $r=zeroes(3,3);
	my $cos=$angle->cos;
	my $sin=$angle->sin;
	$r->diagonal(0,1).=$cos + ($u * $u) * (1-$cos);
	$r(1,0).=$u(0)*$u(1)*(1-$cos)-$u(2)*$sin;
	$r(2,0).=$u(0)*$u(2)*(1-$cos)+$u(1)*$sin;
	$r(0,1).=$u(1)*$u(0)*(1-$cos)+$u(2)*$sin;
	$r(2,1).=$u(1)*$u(2)*(1-$cos)-$u(0)*$sin;
	$r(0,2).=$u(2)*$u(0)*(1-$cos)-$u(1)*$sin;
	$r(1,2).=$u(2)*$u(1)*(1-$cos)+$u(0)*$sin;
	return($r);
}

sub rotation_matrix_vec_vec {
	my ($v1,$v2)=@_;
	my $angle=acos(inner($v1,$v2));
	return(identity(3)) if (abs($angle) < 0.000001);
	my $axis;
	if (abs(abs($angle)-$PI) < 0.01) {
		# avoid numerical instability in case of 180deg rotation
		$axis=zeroes(3);
		my ($i_nz,$i_z)=which_both(abs($v1) > 0.00001);
		if ($i_z->dim(0) == 2) {
			$axis($i_z)->(0).=1;
		} elsif ($i_z->dim(0) == 1) {
			# norm is kept, as v1 is normalized already
			$axis($i_nz).=pdl(1,-1)*$v1($i_nz)->rotate(1);
		} else {
			$axis(pdl(1,2)).=pdl(1,-1)*$v1(pdl(1,2))->rotate(1);
			# renormalize, as we are 'missing' the $v1(0) component
			$axis/=sqrt(inner($axis,$axis));
		}
	} else {
		my $r1_cp=crossp($v1,$v2);
		$axis=$r1_cp/sin($angle);
	}
	my $an=inner($axis,$axis);
	my $r=rotation_matrix_axis_angle($axis,-$angle);
	return($r);
}

sub rotation_matrix_vec_vec_axis {
	my ($v1,$v2,$axis)=@_;
	my $v1a = inner($v1,$axis);
	my $v2a = inner($v2,$axis);
	if (abs($v1a-$v2a) > 0.0001) {
		warn "cannot rotate, vectors are not in the same plane";
		return(undef);
	}
	my $v1p = $v1-($v1a*$axis);
	my $v2p = $v2-($v2a*$axis);
	$v1p/=sqrt(inner($v1p,$v1p));
	$v2p/=sqrt(inner($v2p,$v2p));
	my $angle=acos(inner($v1p,$v2p));
	return(identity(3)) if ($angle < 0.000001);
	my $cp=crossp($v1p,$v2p);
	if (inner($axis,$cp) > 0) {
		$angle=-$angle;
	}

	return(rotation_matrix_axis_angle($axis,$angle));
}

our (@ISA, @EXPORT_OK);
BEGIN {
	require Exporter;
	@ISA = qw(Exporter);
	@EXPORT_OK = qw(&rotation_matrix_axis_angle &rotation_matrix_vec_vec &rotation_matrix_vec_vec_axis);  # symbols to export on request
}
1;
