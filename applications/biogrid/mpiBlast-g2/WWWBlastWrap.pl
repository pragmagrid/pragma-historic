#!/usr/bin/perl

# Global Variables
$scratch_space = "/netapp/blastdb/scratch";
$MPIBLASTCONF = "/usr/local/etc/mpiblast.conf";

# Run as pbsuser
#$> = 501;

#if($> != 501){
#	HtmlError(500, "server Error", "Unable to setuid to pbsuser");
#}

# Check to see if the lock file exists (which means that the 
# blast databases are being updated). Note that the name of the
# lock file must be the same as the name of the lock file in
# the script that updates the databases (currently /etc/cron.weekly/updateblastdb.pl).
$lockfile = '/netapp/blastdb/scratch/blastdb.lock';

if(-e $lockfile){
	HtmlError(500, "Please come back later", "Blast databases are being updated");
}

# Provide a CGI wrapper to the mpiBLAST program
%data = ParseFormData();

# Provide additional keys
$data{'JOB_ID'} = $ENV{'REMOTE_ADDR'};

$data{'TMP_SEQ_FILE'} = CreateUniqueFilename("$scratch_space/$data{'JOB_ID'}.", ".seq");
$data{'TMP_OUTPUT_FILE'} = CreateUniqueFilename("$scratch_space/$data{'JOB_ID'}.", ".html");
$data{'ERROR_LOG_FILE'} = CreateUniqueFilename("$scratch_space/$data{'JOB_ID'}.", ".error");
$data{'STDOUT_LOG_FILE'} = CreateUniqueFilename("$scratch_space/$data{'JOB_ID'}.", ".stdout");

# Set the number of processors to use. There are hardcodes because the databases must be 
# preformatted before the user can perform a blast search.
if($data{'DATALIB'} eq "nt"){
	$data{'NUMPROC'} = 20;
}

if($data{'DATALIB'} eq "nr"){
	$data{'NUMPROC'} = 20;
}

if($data{'DATALIB'} eq "pdb"){
	$data{'NUMPROC'} = 10;
}

# Debug form input
#$@foo = keys %data;
# 
#HtmlHeader("Debug");
#foreach $key (@foo){
#	print "$key=$data{$key}<br>\n";
#}
#HtmlFooter();


# Validate the users input
ValidateFormData(\%data);

# Save the sequence to a temp file
SaveSeq(\%data);

# Run mpiBLAST
Blast(\%data);

# Keep looping until the blast job is done
while(-e $data{'TMP_SEQ_FILE'}){
	sleep 10;
}

# The blast job has finished!
if(!HtmlResults(\%data)){
	DumpErrorLogToHtml($data{'ERROR_LOG_FILE'});

	if(-e $data{'ERROR_LOG_FILE'}){
        	unlink($data{'ERROR_LOG_FILE'});
	}
 
	if(-e $data{'STDOUT_LOG_FILE'}){
        	unlink($data{'STDOUT_LOG_FILE'});
	}

}

exit;

sub ParseFormData()
{
	my %form;
	
	my ($query_string, $requestMethod, @key_value_pairs, $key_value, @tmp_keys);

	$requestMethod = $ENV{'REQUEST_METHOD'};
	
	if(!defined($requestMethod)){
		HtmlError(500, "Server Error", "No request method specified");
	}

	if($requestMethod eq "GET"){
		$query_string = $ENV{'QUERY_STRING'};
	}
	else{
		if ($requestMethod eq "POST"){
			read(STDIN, $query_string, $ENV{'CONTENT_LENGTH'});
		}
		else{
			HtmlError(500, "Server Error", "Server uses unsupported method");
		}
	}

	@key_value_pairs = split(/Content-Disposition: form-data; name=/, $query_string);

	$query_string =~ s/\n/<br>/g;

	foreach $key_value (@key_value_pairs) {

		if($key_value =~ /"(\w+)"\s+(.*)\s+-+/s){
			$form{$1} = $2;
		}
		else{ # Grab the file name
			if($key_value =~ /"(\w+)";\s+(.*)\s+-+/s){
                        	$form{$1} = $2;
                	}
		}
	}

	@tmp_keys = keys %form;
 
	foreach $key (@tmp_keys){
		# We'll add our own end-of-line characters
        	chop($form{$key});
			
		# Make sure we work with windows
		$form{$key} =~ s/\r//sg; 
	}

	return %form;
}

sub ValidateFormData(\%)
{
	my $data_ref = shift;

	# Make sure that the program, query and database are all consistent
	if($data_ref->{'PROGRAM'} eq "blastp"){
		# Must be applied to a protein database
		if($data_ref->{'DATALIB'} eq "nt"){
			HtmlError(500, "Input Error", 
				"blastp must be used with a database of amino acid sequences");
		}
	}

	if($data_ref->{'PROGRAM'} eq "blastn"){ 
		# Must be applied to a nucleotide database
                if($data_ref->{'DATALIB'} ne "nt"){
                        HtmlError(500, "Input Error",  
                                "blastn must be used with a database of nucleotide sequences");
                }
        }

	if($data_ref->{'PROGRAM'} eq "blastx"){ 
		# Must be applied to a protein database
                if($data_ref->{'DATALIB'} eq "nt"){
                        HtmlError(500, "Input Error",  
                                "blastx must be used with a database of amino acid sequences");
                }
        }

	if($data_ref->{'PROGRAM'} eq "tblastn"){ 
		# Must be applied to a nucleotide database
                if($data_ref->{'DATALIB'} ne "nt"){
                        HtmlError(500, "Input Error",
                                "tblastn must be used with a database of nucleotide sequences");
                }
        }

	if($data_ref->{'PROGRAM'} eq "tblastx"){ 
		# Must be applied to a nucleotide database
                if($data_ref->{'DATALIB'} ne "nt"){ 
                        HtmlError(500, "Input Error",
                                "tblastx must be used with a database of nucleotide sequences");
                }
        }
}

sub HtmlError($$$)
{
	my $error_code = shift;
	my $error_header = shift;
	my $error_msg = shift;

	#print "Content-type: text/html\n";
	#print "Status: $error_code, $error_header\n";

	print <<End_of_Error;

Content-type: text/html

<HTML>
<HEAD>
	<TITLE>BLAST - Unexpected Error</TITLE>
</HEAD>
<BODY text="#000000" bgcolor="#FFFFFF" link="#0000EE" vlink="#551A8B" alink="#FF0000">
<H1>$error_header</H1>
<HR>$error_msg</HR>
</BODY>
</HTML>

End_of_Error

	exit(1);	
}

sub HtmlHeader($)
{
	my $title = shift;

	print "Content-type: text/html\n\n";	
	print "<HTML>\n<HEAD>\n\t<TITLE>$title</TITLE>\n</HEAD>\n" .
	      "<BODY text="#000000" bgcolor="#FFFFFF" link="#0000EE" vlink="#551A8B" alink="#FF0000">\n";
}

sub HtmlFooter()
{
	print "</BODY>\n</HTML>\n";
}

sub BlastCommandLine(\%)
{
	# Construct the Blast command line from the
	# user specifed options.
	my $data_ref = shift;
	
	$blast_matrix = $data_ref->{'MAT_PARAM'};
	($blast_matrix) = split /\s+/, $blast_matrix, 0; 
	
	if($data_ref->{'FILTER'} eq "L"){
		$blast_filter = "T";
	}
	else{
		if($data_ref->{'FILTER'} eq "m"){
			$blast_filter = "m S";
		}
		else{
			$blast_filter = "F";
		}
	}
	
	if(defined($data_ref->{'UNGAPPED_ALIGNMENT'})){
		$blast_gapped_alignment = "F";
	}
	else{
		$blast_gapped_alignment = "T";
	}
	
	if($data_ref->{'GENETIC_CODE'} =~ /\((\d)\)/){
		$blast_genetic_code = $1;	
	}
	else{
		$blast_genetic_code = 1;
	}
	
	$blast_advance_opt = 
		ParseAdvancedOpt($data_ref->{'PROGRAM'}, $data_ref->{'OTHER_ADVANCED'});

	# Create the command line to pass to mpiBlast
	my $c_line = "-d $data_ref->{'DATALIB'}$data_ref->{'NUMPROC'} " .
		     "-p $data_ref->{'PROGRAM'} " .
		     "-i $data_ref->{'TMP_SEQ_FILE'} " .
		     "-o $data_ref->{'TMP_OUTPUT_FILE'} " .
		     "-T T " .
		     "-M $blast_matrix " .
		     "-F $blast_filter " . 
		     "-e $data_ref->{'EXPECT'} " .
		     "-w $data_ref->{'OOF_ALIGN'} " .
		     "-g $blast_gapped_alignment " .
		     "-v $data_ref->{'DESCRIPTIONS'} " .
		     "-b $data_ref->{'ALIGNMENTS'} ". 
		     "-m $data_ref->{'ALIGNMENT_VIEW'} " .
		     "-D $blast_genetic_code " .
		     "$blast_advance_opt ";

	return $c_line;
}

sub ParseAdvancedOpt($$)
{
	# Parse the user supplied command line options
	my $blast_prog = shift;
	my $opts = shift;

	my (%opt_hash, @key_pair, $key, $value, $c_line);

	$c_line = "";
	
	@key_pair = split /\s+/, $opts;

	if((scalar(@key_pair) % 2) != 0){
		HtmlError(500, "Input Error", "Advanced options are missing an entry: $opts");
	}

	while(scalar(@key_pair) > 0){
		$value = pop @key_pair;
		$key = pop @key_pair;
		
		# Make sure all values at least look like numbers!
		if(($value =~ /^[-|\+]?\d*.?\d*$/) && !($value =~ /[a-zA-Z]/)){
			 $opt_hash{$key} = $value;
		}
		else{
			HtmlError(500, "Input Error", "Error reading numeric value from " .
				"\"Other advanced options\": $opts");
		}
	}

	# Note that the @key_pair and $key variables are reused in the
	# next bit of code for no good reason.
	if($blast_prog eq "blastn"){
		@key_pair = ("-G", "-E", "-q", "-r", "-W");
		
		foreach $key (@key_pair){
			if(defined($opt_hash{$key})){
				$c_line .= "$key $opt_hash{$key} ";
			}
		}	
	}
	else{ # blastp, blastx, tblastn
		@key_pair = ("-G", "-E", "-W");
 
                foreach $key (@key_pair){
                        if(defined($opt_hash{$key})){
                                $c_line .= "$key $opt_hash{$key} ";
                        }
                }
	}

	return $c_line;
}

sub Blast(\%)
{
 	my $data_ref = shift;
 
 	# Note that we must allocate nodes for NUMPROC + 1 workers (a worker for each of the
	# NUMPROC + 1 fragments. This is a counting peculiairty of mpiBLAST) + 1 node for the master
	# process. Hence a total of NUMPROC + 2 processors are required.
        $NUMNODES = $data_ref->{'NUMPROC'} + 2;
        $PBSINFO = $data_ref->{'JOB_ID'};
	
	my $c_line = BlastCommandLine(%$data_ref);
	
	# Create a temporary file that we'll pass to OpenPBS via qsub
	$tmpFile = CreateUniqueFilename("$scratch_space/blast_script.", ".pl");
 
	open (SCRIPTFILE, ">$tmpFile") or 
		HtmlError(500, "Blast Error", "Unable to create temporary file: $tmpFile");
 
	# Write the contents of the script file
	# Note that the PBS -N command is limited to a string of no more than
	# fifteen characters.

	print SCRIPTFILE '#!/usr/bin/perl'."\n\n";
	print SCRIPTFILE '#PBS -N mpiBLAST'."\n";
	print SCRIPTFILE '#PBS -o '."$data_ref->{'STDOUT_LOG_FILE'}\n";
	print SCRIPTFILE '#PBS -e '."$data_ref->{'ERROR_LOG_FILE'}\n\n";
	print SCRIPTFILE 'if(-e $ENV{PBS_NODEFILE} ){'."\n";
	print SCRIPTFILE '      # Count the number of nodes available to this'."\n";
	print SCRIPTFILE '      # process.'."\n";
	print SCRIPTFILE '      $np = CountNodes($ENV{PBS_NODEFILE});'."\n";
	print SCRIPTFILE '      if($np < 1 || $np > '.$NUMNODES.'){'."\n";
	print SCRIPTFILE '              die "Error: Could not determine number of available nodes.\n";'."\n";
	print SCRIPTFILE '      }'."\n";
	print SCRIPTFILE '      # Make sure we recieved the requested number of nodes'."\n";
	print SCRIPTFILE '      if($np < '.$NUMNODES.'){'."\n";
	print SCRIPTFILE '              die "Error: Did not receive the requested number of nodes\n".'."\n";
	print SCRIPTFILE '                  "Asked for '.$NUMNODES.', received $np nodes.\n";'."\n";
	print SCRIPTFILE '      }'."\n";
	print SCRIPTFILE "\n";
	print SCRIPTFILE '      # Fire up mpiBLAST'."\n";
	print SCRIPTFILE '      system("lamboot $ENV{PBS_NODEFILE} > /dev/null");'."\n";
	print SCRIPTFILE '      $wallClockTime = time;'."\n";
	print SCRIPTFILE '      system("mpirun -np $np mpiblast --removedb".'."\n";
	print SCRIPTFILE '             " --config-file='.$MPIBLASTCONF.' '.$c_line.'");'."\n"; 
	print SCRIPTFILE '      $wallClockTime = time - $wallClockTime;'."\n";
	print SCRIPTFILE '      unlink("'.$data_ref->{'TMP_SEQ_FILE'}.'");'."\n";
	print SCRIPTFILE '      system("lamclean > /dev/null");'."\n";
	print SCRIPTFILE '      system("lamhalt > /dev/null");'."\n";
	print SCRIPTFILE '      print "Wall clock time used is " . $wallClockTime . " sec\n";'."\n";
	print SCRIPTFILE "\n";
	print SCRIPTFILE "}\n";
	print SCRIPTFILE 'else{'."\n";
	print SCRIPTFILE '      unlink("'.$data_ref->{'TMP_SEQ_FILE'}.'");'."\n";
	print SCRIPTFILE '      die "Error: Unable to find the nodefile\n";'."\n";
	print SCRIPTFILE "}\n";
	print SCRIPTFILE "\n";
	print SCRIPTFILE "\n";
	print SCRIPTFILE 'sub CountNodes($)'."\n";
	print SCRIPTFILE "{\n";
	print SCRIPTFILE '      my $filename = shift;'."\n";
	print SCRIPTFILE "\n";
	print SCRIPTFILE '      my ($line, $np);'."\n";
	print SCRIPTFILE "\n";
	print SCRIPTFILE '      $np = 0;'."\n";
	print SCRIPTFILE "\n";
	print SCRIPTFILE '      open(NODEFILE, $filename) or die "Error: Unable to open $filename\n";'."\n";
	print SCRIPTFILE "\n";
	print SCRIPTFILE '      LINE: while($line = <NODEFILE>){'."\n";
	print SCRIPTFILE '              # Just to be safe, ignore any lines that start with'."\n";
	print SCRIPTFILE '              # a comment character (i.e. #).'."\n";
	print SCRIPTFILE "\n";
	print SCRIPTFILE '              next LINE if $line =~ /^#/;'."\n";
	print SCRIPTFILE "\n";
	print SCRIPTFILE '              $np ++;'."\n";
	print SCRIPTFILE '      }'."\n";
	print SCRIPTFILE "\n";
	print SCRIPTFILE '      close(NODEFILE);'."\n";
	print SCRIPTFILE "\n";
	print SCRIPTFILE '      return $np;'."\n";
	print SCRIPTFILE "}\n";
	print SCRIPTFILE "\n";
	 
	close(SCRIPTFILE);

	# Submit our job to OpenPBS
	if(system("/usr/local/bin/qsub -l nodes=$NUMNODES $tmpFile > /dev/null")){
	        # qsub returns > 0 on an error. Clean up the temporary files
		unlink($tmpFile);
		unlink($data_ref->{'TMP_SEQ_FILE'});
	        HtmlError(500, "Blast Error", "Error executing qsub.<br>".
	               "Command was: qsub -l nodes=$NUMNODES $tmpFile");
	}

	unlink($tmpFile);
 
	# All done!
}
	 
sub CreateUniqueFilename($$)
{
        my $basename = shift;
	my $extension = shift;
	my $i = 0;

        while( -e ($basename . $i . $extension) ){
                $i++;
        }
 
        return $basename . $i . $extension;
}

sub SaveSeq(\%)
{
	my $data_ref = shift;
	my $seq_file;


	$seq_file = $data_ref->{'TMP_SEQ_FILE'};

	# Open the seq file for wirting
	open(SEQFILE, ">$seq_file") or 
		HtmlError(500, "IO Error", "Unable to open $seq_file to write seq data to disk");

	# Is the data in SEQUENCE or SEQFILE?
	if(length($data_ref->{'SEQUENCE'}) > 0){
		# Write the sequence to disk
		print SEQFILE $data_ref->{'SEQUENCE'} . "\n";
	}
	else{
		# Windows Explorer clients will append a "Content-Type" string that
		# we must remove.
		if($data_ref->{'SEQFILE'} =~ /Content-Type:\s\S+\s+(.*)/s){
			$data_ref->{'SEQFILE'} = $1;
		}

		if(length($data_ref->{'SEQFILE'}) > 0){
			# Write seq to disk
			print SEQFILE $data_ref->{'SEQFILE'} . "\n";
		}
		else{
			unlink($seq_file);
			HtmlError(500, "Input Error", "No input sequence found!");
		}
	}	
	
	close(SEQFILE);
}

sub HtmlResults(\%)
{
	my $data_ref = shift;
	my ($line, @path, $results_file);

	# Check the size of the output file
	if(-z $data_ref->{'TMP_OUTPUT_FILE'}){
		return 0;
	}
	
	@path = split /\//, $data_ref->{'TMP_OUTPUT_FILE'};
        $results_file = pop @path;

	if(symlink($data_ref->{'TMP_OUTPUT_FILE'}, "/home/apache/html/blast/BlastResults/$results_file") == 0){
		 HtmlError(500, "IO error", "Unable to copy blast results!");
	}

	# Note that there is a problem cleaning up the output files from PBS.
	# By the time we get here (in this script) the stdout and stderror files
	# may not be output yet.
	
	# Point the user to this semi-persistent file
	print "Location: https://jojo.lanl.gov/blast/BlastResults/$results_file\n\n";
	
	return 1;

	# The code below (which never gets executed) is the old way of returning the
	# results. The dyanmic web page was too easy to lose.

	# Open the mpiBlast output file and write the contents to STDOUT
	open(INPUT, $data_ref->{'TMP_OUTPUT_FILE'}) or
		return 0;

	while($line = <INPUT>){
		print $line;
	}

	close(INPUT);

	return 1;
}

sub DumpErrorLogToHtml($)
{
	my $error_file = shift;
	my $line;

	# Make sure that the error file exists
	if(!(-e $error_file)){
		HtmlError(500, "IO Error", "Vauge and unsettling error \# 1 has occured. Blame Jason.");
	}

	# Does the error log contain anything?
	if(-z $error_file){
		HtmlError(500, "IO Error", "Vauge and unsettling error \# 2 has occured. Blame Jason.");
	}

	open(INPUT, $error_file) or
		HtmlError(500, "IO Error", "Unable to open error log");

	HtmlHeader("Blast Error Log");

	while($line = <INTPUT>){
		$line =~ tr/\n/<br>/;

		print STDOUT $line;
	}

	HtmlFooter();

	close(INPUT);
} 
