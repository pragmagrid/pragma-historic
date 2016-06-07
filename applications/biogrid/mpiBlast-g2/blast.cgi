#!/bin/csh -f

#
# $Id: blast.cgi,v 1.1.1.2 2005/01/27 03:02:20 cwwang Exp $
#

#echo "Content-type: text/html"
#echo ""

#setenv DEBUG_COMMAND_LINE TRUE
setenv BLASTDB db

#./blast.REAL
./WWWBlastWrap.pl
