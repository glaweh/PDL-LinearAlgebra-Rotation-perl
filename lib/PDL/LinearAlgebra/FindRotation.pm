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
	my %opts=%{pop()} if (ref $_[-1] eq 'HASH');
	my ($c1,$c2,@constraints)=@_;
	my $eps_length=0.0001;
	my $eps_vol=0.1;
	my $eps_angle=0.01;
	$eps_length=$opts{eps_length} if (defined $opts{eps_length});
	$eps_vol   =$opts{eps_vol}    if (defined $opts{eps_vol});
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
		print STDERR "c1v: " . det($c1) . "\n";
		print STDERR "c2v: " . det($c2) . "\n";
	}

	my $c1va = $c1v->abs;
	my $c2va = $c2v->abs;
	if (abs($c1va-$c2va)>$eps_vol) {
		print STDERR "find_rotation: conv_cell volume mismatch: $c1va / $c2va\n";
		return(undef);
	} else {
		print STDERR "find_rotation: conv_cell volume match: $c1va\n" if ($DEBUG > 0);
	}

	my $c12_match=long(abs($c1l-$c2l->dummy(0)) < $eps_length);
	print_matrix('c12_match',$c12_match) if ($DEBUG > 1);

	my $c12_free=$c12_match->copy;
	foreach my $constraint (@constraints) {
		if ($c12_match(@{$constraint}) != 1) {
			print STDERR "Failed fulfilling constraint: c1l($constraint->[0]) == c2l($constraint->[1])\n";
			print STDERR "  " . $c1l($constraint->[0]) . " != " . $c2l($constraint->[1]) . "\n";
			return(undef);
		} elsif ($DEBUG > 2) {
			print STDERR "Fulfilled constraint: c1l($constraint->[0]) == c2l($constraint->[1])\n";
		}
		$c12_free($constraint->[0],:).=0;
		$c12_free(:,$constraint->[1]).=0;
		$c12_free(@{$constraint}).=1;
	}
	print_matrix('c12_free',$c12_free) if ($DEBUG > 1);
	my @i_0=which($c12_free(0,:;-))->list;
	my @i_1=which($c12_free(1,:;-))->list;
	my @try_axis_map;
	foreach my $i_0 (@i_0) {
		foreach my $i_1 (@i_1) {
			next if ($i_0 == $i_1);
			push @try_axis_map,[ 0, $i_0, 1, $i_1 ];
		}
	}
	if ($#try_axis_map < 0) {
		print STDERR "No possible axis mapping found!\n";
		return(undef);
	}

	my @rotation;
	foreach my $axis_map (@try_axis_map) {
		print STDERR "trying axis mapping $axis_map->[0]:$axis_map->[1] $axis_map->[2]:$axis_map->[3]\n" if ($DEBUG > 1);
		foreach my $hand_prefactor (1,-1) {
			print STDERR "  handedness prefactor: $hand_prefactor\n" if ($DEBUG > 1);
			my $rotation = rotation_sequence($hand_prefactor*$c1,$c2,@{$axis_map});
			if (defined $rotation) {
				$rotation*=$hand_prefactor;
				my $c1r = $c1 x $rotation;
				if (check_vector_match('    find_rotation',$c1r,$c2)) {
					print_matrix("    Found rotation: ",$rotation,'    ') if ($DEBUG > 0);
					push @rotation,$rotation;
				} elsif ($DEBUG > 0) {
					print_matrix("    Invalid rotation: ",$rotation,'    ') if ($DEBUG > 0);
					print_matrix("      rotated: ",$c1r,'      ');
					print_matrix("       target: ",$c2,'      ');
				}
			} elsif ($DEBUG > 1) {
				print STDERR "    No rotation found.\n";
			}
		}
	}
	my $nrot = $#rotation +1;
	unless ($nrot > 0) {
		print STDERR "find_rotation: not found\n" if ($DEBUG > 0);
		return(undef);
	}
	print_matrix("find_rotation: success, found $nrot possible rotations!",cat(@rotation),'  ') if ($DEBUG > 0);
	return(\@rotation);
}

our (@ISA, @EXPORT_OK);
BEGIN {
	require Exporter;
	@ISA = qw(Exporter);
	@EXPORT_OK = qw(&find_rotation);  # symbols to export on request
}
1;
