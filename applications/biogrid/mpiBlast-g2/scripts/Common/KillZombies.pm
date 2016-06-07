#!/usr/bin/perl
package Common::KillZombies;
require Exporter;

our @ISA    = qw(Exporter);
our @EXPORT = qw(killZombies
                 getRemoteTemp
                 cleanUpTemp);         

use IPC::Open3;
use FileHandle;

sub killZombies {

    my $bhost = shift;

    my $exec = "ssh $bhost \'killall -9 ".$ENV{'BLASTG2_LOCATION'}."\/bin\/mpiblast-g2\'";

    print $exec."\n";

    ## Kill the zombie job 
    my ($wrt,$rdr,$err) = (FileHandle->new,FileHandle->new,FileHandle->new);

    my $pid = open3($wrt,$rdr,$err,$exec);  # Execution the command

    if($pid) {

        waitpid($pid,0);

        while( $line=<$rdr> ) {

            chomp($line);

            $outMsg.=$line;
        }

        while( $line=<$err> ) {
            $errMsg.=$line;
        }
    }

    close($wrt);
    close($rdr);
    close($err);

}

sub getRemoteTemp {

    my $bhost = shift;

    my $exec = "ssh $bhost \'cat -b ".$ENV{'BLASTG2_LOCATION'}."\/etc\/mpiblast-gasscopy.conf \'";

    print $exec."\n";

    my $tmp;

    my ($wrt,$rdr,$err) = (FileHandle->new,FileHandle->new,FileHandle->new);

    my $pid = open3($wrt,$rdr,$err,$exec);  # Execution the command

    if($pid) {

        waitpid($pid,0);

        my $cnt = 0;

        while( $line=<$rdr> ) {

            chomp($line);

            if( $line =~ m/^\s+(3)\s+(.*)/) { $tmp = $2; }

            $outMsg.=$line;
        }

        while( $line=<$err> ) {
            $errMsg.=$line;
        }
    }

    close($wrt);
    close($rdr);
    close($err);

    return $tmp;
}

sub cleanUpTemp {

    my ($bhost,$tmp) = @_;

    my $exec = "ssh $bhost \" ls -l $tmp | grep $ENV{'USER'} | awk '{print \\\"$tmp/\\\" \\\$9}' | xargs rm -rf \"";

    print $exec."\n";

    my ($wrt,$rdr,$err) = (FileHandle->new,FileHandle->new,FileHandle->new);

    my $pid = open3($wrt,$rdr,$err,$exec);  # Execution the command

    if($pid) {

        waitpid($pid,0);

        while( $line=<$rdr> ) {
            $outMsg.=$line;
        }

        while( $line=<$err> ) {
            $errMsg.=$line;
        }
    }

    close($wrt);
    close($err);
    close($rdr);

    return 1;

}

1;
