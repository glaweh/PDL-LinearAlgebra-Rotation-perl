package PDL::LinearAlgebra::FindRotation;
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::LinearAlgebra::Rotation qw{rotation_matrix_vec_vec rotation_matrix_vec_vec_axis};

our $DEBUG=0;
my $PI=2*acos(0);

sub print_matrix {
	my ($prefix,$matrix,$indent)=@_;
	my $m=$matrix->copy;
	$indent='' unless (defined $indent);
	$m->flat->where(abs($m) < 0.00001).=0;
	my $mp = '' . $m . '';
	$mp =~ s/(^|\n)/$1$indent/sg;
	print STDERR $prefix . $mp . "\n";
}

sub check_vector_match {
	my ($prefix,$c1,$c2)=@_;
	my $match=whichND(($c1->dummy(1)-$c2)->abs->sumover < 0.001);
	my $c2_to_c1=zeroes(long,3)-1;
	$c2_to_c1($match(0,:)).=$match(1,:);
	print STDERR "$prefix: c2_to_c1: $c2_to_c1\n" if ($DEBUG > 1);
	if (sum($c2_to_c1 >= 0) < 3) {
		print STDERR "$prefix: rotated cell does not match target values:\n  c2_to_c1: $c2_to_c1\n" if ($DEBUG > 0);
		return(0);
	} else {
		print STDERR "$prefix: rotated cell matches target values\n" if ($DEBUG > 0);
		return(1);
	}
}

sub rotation_sequence {
	my ($c1,$c2,@axis)=@_;
	# match $c1(:,$axis[0]), $c2(:,$axis[1])
	# rotate around $c2(:,$axis[1])
	#   to match $c1(:,$axis[2]), $c2(:,$axis[3])
	my ($c1l,$c1n,$c1a);         # auxiliary c1 data
	my ($c2l,$c2n,$c2a);         # auxiliary c2 data
	my $c1_orig=$c1;             # c1 original matrix
	my $rotation=identity(3);    # accumulator for the rotation matrix
	$c1l = sqrt(inner($c1,$c1));
	$c2l = sqrt(inner($c2,$c2));
	$c1n = $c1/$c1l->dummy(0);
	$c2n = $c2/$c2l->dummy(0);
	$c1a = acos(inner($c1n(:,pdl(long,1,2,0)),$c1n(:,pdl(long,2,0,1))))*180/$PI;
	$c2a = acos(inner($c2n(:,pdl(long,1,2,0)),$c2n(:,pdl(long,2,0,1))))*180/$PI;

	# get the first rotation
	my $r=rotation_matrix_vec_vec($c1n(:,$axis[0];-),$c2n(:,$axis[1];-));

	print_matrix("free rotation ($axis[0],$axis[1]): ",$r) if ($DEBUG > 1);
	$rotation = $rotation x $r;
	$c1 = $c1_orig x $rotation;
	$c1n = $c1/$c1l->dummy(0);
	print_matrix("c1_rotated: ",$c1) if ($DEBUG > 1);

	$r=rotation_matrix_vec_vec_axis($c1n(:,$axis[2];-),$c2n(:,$axis[3];-),$c1n(:,$axis[0];-));
	unless (defined $r) {
		return(undef);
	}
	print_matrix("rotation around $axis[0] ($axis[2],$axis[3]): ",$r) if ($DEBUG > 1);
	$rotation = $rotation x $r;
	$c1 = $c1_orig x $rotation;
	$c1n = $c1/$c1l->dummy(0);
	print_matrix("c1_rotated: ",$c1) if ($DEBUG > 1);

	return(undef) unless check_vector_match('rotation_sequence',$c1,$c2);
	return($rotation);
}

=head2 find_rotation

Given two sets of 3 vectors each, where one set is a (possibly improperly)
rotated version of the other one, find_rotation() finds the corresponding
rotation matrix.

In each of the two parallelepipeds, one axis and/or one face may be
flagged as "special", which is especially useful if there are vectors with
the same length, but inequivalent under other aspects (think crystallographic
axis).

=cut

sub find_rotation {
	# finds the rotation between two conventional cells
	# may contain a flip to match handedness
	# result: matrix M: c2 = M c1
	my ($c1,$c1s,$c1c,$c2,$c2s,$c2c,$opts)=@_;
	my $eps_length=0.0001;
	my $eps_vol=0.1;
	my $eps_angle=0.01;
	if (defined $opts and ref $opts eq 'HASH') {
		$eps_length=$opts->{eps_length} if (defined $opts->{eps_length});
		$eps_vol   =$opts->{eps_vol}    if (defined $opts->{eps_vol});
	}
	my ($c1l,$c1n,$c1a,$c1v);
	my ($c2l,$c2n,$c2a,$c2v);
	$c1l = sqrt(inner($c1,$c1));
	$c2l = sqrt(inner($c2,$c2));
	$c1n = $c1/$c1l->dummy(0);
	$c2n = $c2/$c2l->dummy(0);
	$c1a = acos(inner($c1n(:,pdl(long,1,2,0)),$c1n(:,pdl(long,2,0,1))))*180/$PI;
	$c2a = acos(inner($c2n(:,pdl(long,1,2,0)),$c2n(:,pdl(long,2,0,1))))*180/$PI;
	$c1v = inner($c1(:,0;-),crossp($c1(:,1;-),$c1(:,2;-)));
	$c2v = inner($c2(:,0;-),crossp($c2(:,1;-),$c2(:,2;-)));
	if ($DEBUG > 0) {
		print STDERR "c1: $c1\n";
		print STDERR "c2: $c2\n";
		print STDERR "c1l: $c1l\n";
		print STDERR "c2l: $c2l\n";
		print STDERR "c1a: $c1a\n";
		print STDERR "c2a: $c2a\n";
	}

	my $c1va = $c1v->abs;
	my $c2va = $c2v->abs;
	if (abs($c1va-$c2va)>$eps_vol) {
		print STDERR "find_rotation: conv_cell volume mismatch: $c1va / $c2va\n";
		return(undef);
	} else {
		print STDERR "find_rotation: conv_cell volume match: $c1va\n" if ($DEBUG > 0);
	}

	my ($centered,$special_axis);
	# check centering
	if (($centered=(defined $c1c ? 1 : 0)) xor (my $c2=defined $c2c)) {
		print STDERR "find_rotation: conv_cell centering setting mismatch: $centered/$c2\n";
		return(undef);
	}
	if (($special_axis=(defined $c1s ? 1 : 0)) xor (my $s2=defined $c2s)) {
		print STDERR "find_rotation: conv_cell special axis setting mismatch: $special_axis/$s2\n";
		return(undef);
	}
	print STDERR "centered/special: $centered/$special_axis\n" if ($DEBUG > 1);

	my $c1_orig = $c1;
	my ($c1_1,@c2_1_candidates);
	my ($c1_2,@c2_2_candidates);
	if ($special_axis) {
		if (abs($c1l($c1s)-$c2l($c2s)) > $eps_length) {
			print STDERR "find_rotation: special axis have different lengths: " . $c1l($c1s) . " / " . $c2l($c2s) . "\n";
			return(undef);
		}
		$c1_1=$c1s;
		@c2_1_candidates=($c2s);
		if ($centered and ($c1c!=$c1s)) {
			$c1_2 = (grep { ($_ != $c1s) and ($_ != $c1c) } 0 .. 2)[0];
			@c2_2_candidates = grep { ($_ != $c2s) and ($_ != $c2c) } 0 .. 2;
			if (abs($c1l($c1_2)-$c2l($c2_2_candidates[0]))>0.001) {
				print STDERR "find_rotation: axis to-match have different lengths: " . $c1l($c1_2) . " / " . $c2l($c2_2_candidates[0]) . "\n";
				return(undef);
			}
		} else {
			# find by matching lengths
			$c1_2    = (grep { ($_ != $c1s) } 0 .. 2)[0];
			my $c2mask = (abs($c2l-$c1l($c1_2)) < $eps_length);
			$c2mask($c2s).=0;
			@c2_2_candidates = which($c2mask)->list;
		}
	} else {
		my $c2mask;
		$c1_1=0;
		$c2mask = (abs($c2l-$c1l(0)) < $eps_length);
		@c2_1_candidates = which($c2mask)->list;
		$c1_2=1;
		$c2mask = (abs($c2l-$c1l(1)) < $eps_length);
		@c2_2_candidates = which($c2mask)->list;
	}

	if ($#c2_1_candidates < 0) {
		print STDERR "no c2_1_candidates to match\n";
		return(undef);
	} elsif ($DEBUG > 1) {
		print STDERR "c2_1_candidates:\n";
		print STDERR "  $c1_1: " . join(', ',@c2_1_candidates) . "\n";
	}
	if ($#c2_2_candidates < 0) {
		print STDERR "no c2_2_candidates to match\n";
		return(undef);
	} elsif ($DEBUG > 1) {
		print STDERR "c2_2_candidates:\n";
		print STDERR "  $c1_2: " . join(', ',@c2_2_candidates) . "\n";
	}
	my @rotation;
	foreach my $c2_1 (@c2_1_candidates) {
		foreach my $c2_2 (@c2_2_candidates) {
			next if ($c2_1 == $c2_2);
			print STDERR "trying: $c1_1:$c2_1 $c1_2:$c2_2\n" if ($DEBUG > 0);
			my $rotation=rotation_sequence($c1,$c2,$c1_1,$c2_1,$c1_2,$c2_2);
			if (defined $rotation) {
				print_matrix("  Found rotation: ",$rotation,'    ') if ($DEBUG > 0);
				push @rotation,$rotation;
			} elsif ($DEBUG > 1) {
				print STDERR "  No rotation found.\n";
			}

			print STDERR "find_rotation: trying with inverted c1\n" if ($DEBUG > 0);
			$rotation=rotation_sequence(-$c1,$c2,$c1_1,$c2_1,$c1_2,$c2_2);
			if (defined $rotation) {
				$rotation*=-1;
				print_matrix("  Found rotation: ",$rotation,'    ') if ($DEBUG > 0);
				push @rotation,$rotation;
			} elsif ($DEBUG > 1) {
				print STDERR "  No rotation found.\n";
			}
		}
	}
	unless ($#rotation >= 1) {
		print STDERR "find_rotation: not found\n";
		return(undef);
	}
	my $rotation = $rotation[0];
	$c1 = $c1_orig x $rotation;
	if ($DEBUG > 0) {
		print_matrix("rotated: ",$c1);
		print_matrix("target: ",$c2);
		print_matrix("rotation: ",$rotation);
	}
	return(undef) unless check_vector_match('find_rotation',$c1,$c2);
	print STDERR "find_rotation: success!\n" if ($DEBUG > 0);
	return($rotation);
}

our (@ISA, @EXPORT_OK);
BEGIN {
	require Exporter;
	@ISA = qw(Exporter);
	@EXPORT_OK = qw(&find_rotation);  # symbols to export on request
}
1;
