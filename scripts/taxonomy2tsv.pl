#! /usr/bin/env perl

# Original program is from Krona package

# Copyright © 2011, Battelle National Biodefense Institute (BNBI);
# all rights reserved. Authored by: Brian Ondov, Nicholas Bergman, and
# Adam Phillippy
#
# See the LICENSE.txt file included with this software for license information.

use strict;

my ($taxonomy) = @ARGV;
my %data;

# load scientific names for each tax ID

open NAMES, "<$taxonomy/names.dmp" or die "Couldn't open names.dmp";

while ( my $line = <NAMES> )
{
	my ($id, $name, $uniqueName, $class) = split /\t\|\t/, $line;
	
	if ( $class =~ /scientific name/ )
	{
		if ( ! defined $data{$id} )
		{
			$data{$id} = ();
		}
		
		$data{$id}->{'name'} = $name;
	}
}

close NAMES;

# load parents and ranks for each tax ID

open NODES, "<$taxonomy/nodes.dmp" or die "Couldn't open nodes.dmp";

while ( my $line = <NODES> )
{
	$line =~ /(\d+)\t\|\t(\d+)\t\|\t([^\t]+)/;
	
	my $id = $1;
	
	if ( ! defined $data{$id} )
	{
		$data{$id} = ();
	}
	
	$data{$id}->{'parent'} = $2;
	$data{$id}->{'rank'} = $3;
}

close NODES;

foreach my $id ( sort {$a <=> $b} keys %data )
{
	print STDOUT join "\t",
	(
		$id,
		depth($id),
		getParent($id),
		$data{$id}->{'rank'},
		$data{$id}->{'name'}
	);
	print STDOUT "\n";
}

sub depth
{
	my ($id) = @_;
	
	if ( $id == 1 )
	{
		return 0;
	}
	else
	{
		return depth(getParent($id)) + 1;
	}
}

sub getParent
{
	my ($id) = @_;
	return $data{$id}->{'parent'};
}
