use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'PDL::LinearAlgebra::Rotation',
    license             => 'gpl',
    dist_author         => q{Henning Glawe <glawe@web.de>},
    dist_version_from   => 'lib/PDL/LinearAlgebra/Rotation.pm',
    build_requires => {
        'Test::More' => 0,
    },
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
