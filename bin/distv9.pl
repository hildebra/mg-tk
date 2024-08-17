#!/usr/bin/perl -w
use strict;

my $tF = $ARGV[0];
my $tree1 = `cat $tF`;
my @spl = split /\n/,$tree1;# =~ m/(^.*)\n/;
my $tree = $spl[0];
##record the distance of parentheses
my %dis;
my $par = -1;
my @current;
while($tree =~ /./g)
    {if ($& eq '(')
        {$par ++;
        next if $par == 0;
        $current[$#current+1] = $par;
        }
    elsif($& eq ')')
        {(my $tem) = $' =~ /:(\d+\.\d+|\d+)/;
        next if $#current == -1;
        $dis{'node_'.$current[$#current]} = $tem;
        pop @current;
        }
    }

##record the distance of leaves
my @order;
while ($tree =~ /([^\(\):,]+):(\d+\.\d+|\d+)/g)
    {$dis{$1} = $2;
    $order[$#order+1] = $1;
    }

##record parents of leaves
my %pare;
@current = ();
$par = -1;
while($tree =~ /(\(|\)|([^\(\):,]+):)/g)
    {if ($& eq '(')
        {$par ++;
        next if $par == 0;
        $current[$#current+1] = $par;
        }
    elsif($& eq ')')
        {pop @current;
        }
    else{map {$pare{$2}{$_} = 1} @current;
        $pare{$2} = [@current];
        }
    }

##Distance matrix
my %dis2;
foreach my $i (0..$#order)
    {foreach my $j ($i..$#order)
        {if ($i == $j)
            {$dis2{$order[$i]}{$order[$j]} = 0;
            }
        else{my $tem = $dis{$order[$i]} + $dis{$order[$j]};
            my $tem2 = -1;
            foreach my $k (0..$#{$pare{$order[$i]}})
                {last if ($k > $#{$pare{$order[$j]}});
                if ($pare{$order[$i]}[$k] eq $pare{$order[$j]}[$k])
                    {$tem2 = $k;
                    }
                }
            if ($#{$pare{$order[$i]}} != -1)
                {map {$tem += $dis{'node_'.$_}} map {$pare{$order[$i]}[$_]} ($tem2+1)..$#{$pare{$order[$i]}};
                }
            if ($#{$pare{$order[$j]}} != -1)
                {map {$tem += $dis{'node_'.$_}} map {$pare{$order[$j]}[$_]} ($tem2+1)..$#{$pare{$order[$j]}};
                }
            $dis2{$order[$i]}{$order[$j]} = $dis2{$order[$j]}{$order[$i]} = $tem;
            }
        }
    }

##output
print join("\t",'',@order),"\n";
foreach my $i (@order)
    {print join("\t",$i,map {$dis2{$i}{$_}} @order),"\n";
    }
