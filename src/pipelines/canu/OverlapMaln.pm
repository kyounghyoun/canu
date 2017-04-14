
###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' (http://kmer.sourceforge.net)
 #  both originally distributed by Applera Corporation under the GNU General
 #  Public License, version 2.
 #
 #  Canu branched from Celera Assembler at its revision 4587.
 #  Canu branched from the kmer project at its revision 1994.
 #
 #  Modifications by:
 #
 #    Sergey Koren beginning on 2016-FEB-24
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Brian P. Walenz beginning on 2016-MAY-02
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::OverlapMaln;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(malnConfigure malnPrecomputeCheck malnCheck);

use strict;

use File::Path 2.08 qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;
use canu::Gatekeeper;
use canu::HTML;
use canu::Grid_Cloud;

#  Map long reads to long reads with minimap.

sub malnConfigure ($$$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $typ     = shift @_;
    my $bin     = getBinDirectory();

    my $base;                #  e.g., $base/1-overlapper/mhap.sh
    my $path;                #  e.g., $path/mhap.sh

    $base = "correction"  if ($tag eq "cor");
    $base = "trimming"    if ($tag eq "obt");
    $base = "unitigging"  if ($tag eq "utg");

    $path = "$base/1-overlapper";

    caFailure("invalid type '$typ'", undef)  if (($typ ne "partial") && ($typ ne "normal"));

    goto allDone   if (skipStage($asm, "$tag-malnConfigure") == 1);
    goto allDone   if (fileExists("$path/precompute.sh")) && (fileExists("$path/maln.sh"));
    goto allDone   if (fileExists("$path/ovljob.files"));
    goto allDone   if (-e "$base/$asm.ovlStore");
    goto allDone   if (fileExists("$base/$asm.ovlStore.tar"));

    print STDERR "--\n";
    print STDERR "-- OVERLAPPER (maln) (correction)\n"  if ($tag eq "cor");
    print STDERR "-- OVERLAPPER (maln) (trimming)\n"    if ($tag eq "obt");
    print STDERR "-- OVERLAPPER (maln) (assembly)\n"    if ($tag eq "utg");
    print STDERR "--\n";

    make_path($path) if (! -d $path);

    #  Constants.

    my $merSize       = getGlobal("${tag}MalnMerSize");
    my $windowSize       = getGlobal("${tag}MalnWindowSize");

    my $numReads      = getNumberOfReadsInStore($base, $asm);
    my $memorySize    = getGlobal("${tag}MalnMemory");
    my $blockPerGb    = getGlobal("${tag}MalnBlockSize");
    my $blockSize = int($blockPerGb * $memorySize);

    print STDERR "-- Given $memorySize GB, can fit $blockSize reads per block.\n";

    #  Divide the reads into blocks of ovlHashBlockSize.  Each one of these blocks is used as the
    #  table in maln.  Several of these blocks are used as the queries.

    my @blocks;    #  Range of reads to extract for this block
    my @blockBgn;  #  First read in the block
    my @blockLen;  #  Number of reads in the block

    my @hashes;    #  One for each job, the block that is the hash table
    my @convert;   #  One for each job, flags to the maln-ovl conversion program

    push @blocks,   "no zeroth block, makes the loop where this is used easier";
    push @blockBgn, "no zeroth block";
    push @blockLen, "no zeroth block";
    push @hashes,   "no zeroth job";

    for (my $bgn=1; $bgn < $numReads; $bgn += $blockSize) {
        my $end = $bgn + $blockSize - 1;
        $end = $numReads  if ($end > $numReads);

        #print STDERR "BLOCK ", scalar(@blocks), " reads from $bgn through $end\n";

        push @blocks, "-r$bgn-$end";

        push @blockBgn, $bgn;
        push @blockLen, $end - $bgn + 1;
    }

    #  Each maln job will process one block against a set of other blocks.  We'll pick, arbitrarily,
    #  to use num_blocks/4 for that size, unless it is too small.

    my $numBlocks = scalar(@blocks);
    my $qryStride = ($numBlocks < 16) ? (2) : int($numBlocks / 4);

    print STDERR "-- For $numBlocks blocks, set stride to $qryStride blocks.\n";
    print STDERR "-- Logging partitioning to '$path/partitioning.log'.\n";

    open(L, "> $path/partitioning.log") or caExit("can't open '$path/partitioning.log' for writing: $!\n", undef);

    #  Make queries.  Each hask block needs to search against all blocks less than or equal to it.
    #  Each job will search at most $qryStride blocks at once.  So, queries could be:
    #  1:  1 vs 1,2,3  (with self-allowed, and quert block 1 implicitly included)
    #  2:  1 vs 4,5    (with no-self allowed)
    #  3:  2 vs 2,3,4
    #  4:  2 vs 5
    #  5:  3 vs 3,4,5
    #  6:  4 vs 4,5
    #  7:  5 vs 5

    make_path("$path/queries");

    for (my $bid=1; $bid < $numBlocks; $bid++) {

        #  Note that we never do qbgn = bid; the self-self overlap is special cased.

        for (my $qbgn = $bid; $qbgn < $numBlocks; $qbgn += $qryStride) {

            my $qend = $qbgn + $qryStride - 1;                 #  Block bid searches reads in dat files from
            $qend = $numBlocks-1   if ($qend >= $numBlocks);   #  qbgn to qend (inclusive).


            my $job = substr("000000" . scalar(@hashes), -6);  #  Unique ID for this compute

            #  Make a place to save queries.  If this is the last-block-special-case, make a directory,
            #  but don't link in any files.  Without the directory, we'd need even more special case
            #  code down in maln.sh to exclude the -q option for this last block.

            make_path("$path/queries/$job");

            if ($qbgn < $numBlocks) {
                print L "Job ", scalar(@hashes), " computes block $bid vs blocks $qbgn-$qend,\n";

                for (my $qid=$qbgn; $qid <= $qend; $qid++) {
                    my $qry = substr("000000" . $qid, -6);             #  Name for the query block

                    symlink("../../blocks/$qry.fasta", "$path/queries/$job/$qry.fasta");
                }

            } else {
                print L "Job ", scalar(@hashes), " computes block $bid vs itself.\n";
                $qbgn = $bid;  #  Otherwise, the @convert -q value is bogus
            }

            #  This is easy, the ID of the hash.

            push @hashes, substr("000000" . $bid, -6);  #  One new job for block bid with qend-qbgn query files in it
        }
    }

    close(L);

    #  Tar up the queries directory.  Only useful for cloud support.

    runCommandSilently($path, "tar -cf queries.tar queries", 1);
    stashFile("$path/queries.tar");

    #  Create a script to generate precomputed blocks, including extracting the reads from gkpStore.

    #OPTIMIZE
    #OPTIMIZE  Probably a big optimization for cloud assemblies, the block fasta inputs can be
    #OPTIMIZE  computed ahead of time, stashed, and then fetched to do the actual precompute.
    #OPTIMIZE

    open(F, "> $path/precompute.sh") or caFailure("can't open '$path/precompute.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F fetchStoreShellCode("$base/$asm.gkpStore", "$base/1-overlapper", "");
    print F "\n";
    print F getJobIDShellCode();
    print F "\n";
    for (my $ii=1; $ii < scalar(@blocks); $ii++) {
        print F "if [ \$jobid -eq $ii ] ; then\n";
        print F "  rge=\"$blocks[$ii]\"\n";
        print F "  job=\"", substr("000000" . $ii, -6), "\"\n";
        print F "fi\n";
        print F "\n";
    }
    print F "\n";
    print F "if [ x\$job = x ] ; then\n";
    print F "  echo Job partitioning error.  jobid \$jobid is invalid.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "if [ ! -d ./blocks ]; then\n";
    print F "  mkdir -p ./blocks\n";
    print F "fi\n";
    print F "\n";
    print F fileExistsShellCode("./blocks/\$job.fasta");
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "\$bin/gatekeeperDumpFASTQ \\\n";
    print F "  -G ../$asm.gkpStore \\\n";
    print F "  \$rge \\\n";
    print F "  -nolibname \\\n";
    print F "  -noreadname \\\n";
    print F "  -fasta \\\n";
    print F "  -o ./blocks/\$job.input \\\n";
    print F "&& \\\n";
    print F "mv -f ./blocks/\$job.input.fasta ./blocks/\$job.fasta\n";
    print F "\n";

    #  The following lines are ported from OverlapMhap.pm
    print F "if [ ! -e \"$path/blocks/\$job.fasta\" ] ; then\n";
    print F "  echo Failed to extract fasta.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "echo Starting minialign precompute.\n";
    print F "\n";
    print F "#  So minialign writes its output in the correct spot.\n";
    print F "cd $path/blocks\n";
    print F "\n";
    print F "\$bin/minialign \\\n";
    print F "  -NX \\\n";
    print F "  -t ", getGlobal("${tag}malnThreads"), " \\\n";
    print F "  -k $merSize \\\n";
    print F "  -w $windowSize \\\n";
    print F "  -d $path/blocks/\$job.mai.WORKING \\\n";
    print F "  $path/blocks/\$job.fasta \\\n";
    print F "&& \\\n";
    print F "mv -f ./blocks/\$job.mai.WORKING ./blocks/\$job.mai\n";
    print F "\n";
    print F "if [ ! -e \"$path/blocks/\$job.mai\" ] ; then\n";
    print F "  echo minialign failed.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "#  Clean up, remove the fasta input\n";
    print F "rm -f $path/blocks/\$job.fasta\n";
    print F "\n";
    print F "exit 0\n";

    close(F);

    #  Create a script to run maln.

    open(F, "> $path/maln.sh") or caFailure("can't open '$path/maln.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F fetchStoreShellCode("$base/$asm.gkpStore", "$base/1-overlapper", "");
    print F "\n";
    print F getJobIDShellCode();
    print F "\n";
    for (my $ii=1; $ii < scalar(@hashes); $ii++) {
        print F "if [ \$jobid -eq $ii ] ; then\n";
        print F "  blk=\"$hashes[$ii]\"\n";
        print F "  qry=\"", substr("000000" . $ii, -6), "\"\n";
        print F "fi\n";
        print F "\n";
    }

    print F "\n";
    print F "if [ x\$qry = x ]; then\n";
    print F "  echo Error: Job index out of range.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "if [ -e ./results/\$qry.ovb ]; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
 
    print F fetchFileShellCode("$path", "queries.tar", "");
    print F "\n";
    print F "if [ -e ./queries.tar -a ! -d ./queries ] ; then\n";
    print F "  tar -xf ./queries.tar\n";
    print F "fi\n";
    print F "\n";

    print F "if [ ! -d ./results ]; then\n";
    print F "  mkdir -p ./results\n";
    print F "fi\n";
    print F "\n";

    print F fetchFileShellCode("$path", "blocks/\$blk.fasta", "");

    print F "for ii in `ls ./queries/\$qry` ; do\n";
    print F "  echo Fetch blocks/\$ii\n";
    print F    fetchFileShellCode("$path", "blocks/\$ii", "  ");
    print F "done\n";
    print F "\n";

    #  Begin comparison, we loop through query and compare current block to it, if we need to do
    #  self first compare to self, otherwise initialize as empty

    print F "for file in `ls queries/\$qry/*.fasta`; do\n";
    print F "  \$bin/minialign \\\n";
    print F "    -NXxava \\\n";
    print F "    -t ", getGlobal("${tag}malnThreads"), " \\\n";
    print F "    -l ./blocks/\$blk.mai \\\n";
    print F "    \$file \\\n";
    print F "  >> ./results/\$qry.paf.WORKING \n";
    print F "done\n";
    print F "\n";
    print F "mv  ./results/\$qry.paf.WORKING  ./results/\$qry.paf\n";
    print F "\n";
    print F "if [   -e ./results/\$qry.paf -a \\\n";
    print F "     ! -e ./results/\$qry.ovb ] ; then\n";
    print F "  \$bin/mmapConvert \\\n";
    print F "    -G ../$asm.gkpStore \\\n";
    print F "    -o ./results/\$qry.paf.ovb.WORKING \\\n";
    print F "    ./results/\$qry.paf \\\n";
    print F "  && \\\n";
    print F "  mv ./results/\$qry.paf.ovb.WORKING ./results/\$qry.paf.ovb\n";
    print F "fi\n";
    print F "\n";

    if (getGlobal('saveOverlaps') eq "0") {
        print F "if [   -e ./results/\$qry.paf -a \\\n";
        print F "       -e ./results/\$qry.paf.ovb ] ; then\n";
        print F "  rm -f ./results/\$qry.paf\n";
        print F "fi\n";
        print F "\n";
    }

    print F "if [ -e ./results/\$qry.paf.ovb ] ; then\n";
    if (getGlobal("${tag}ReAlign") eq "raw") {
        print F "  \$bin/overlapPair \\\n";
        print F "    -G ../$asm.gkpStore \\\n";
        print F "    -O ./results/\$qry.paf.ovb \\\n";
        print F "    -o ./results/\$qry.ovb \\\n";
        print F "    -partial \\\n"  if ($typ eq "partial");
        print F "    -erate ", getGlobal("corOvlErrorRate"), " \\\n"  if ($tag eq "cor");
        print F "    -erate ", getGlobal("obtOvlErrorRate"), " \\\n"  if ($tag eq "obt");
        print F "    -erate ", getGlobal("utgOvlErrorRate"), " \\\n"  if ($tag eq "utg");
        print F "    -memory " . getGlobal("${tag}malnMemory") . " \\\n";
        print F "    -t " . getGlobal("${tag}malnThreads") . " \n";
    } else {
        print F "  mv -f ./results/\$qry.paf.ovb    ./results/\$qry.ovb\n";
        print F "  mv -f ./results/\$qry.paf.counts ./results/\$qry.counts\n";
    }
    print F "fi\n";

    print F stashFileShellCode("$path", "results/\$qry.ovb",    "");
    print F stashFileShellCode("$path", "results/\$qry.counts", "");
    print F "\n";

    print F "\n";
    print F "exit 0\n";

    close(F);

    if (-e "$path/precompute.sh") {
        my $numJobs = 0;
        open(F, "< $path/precompute.sh") or caFailure("can't open '$path/precompute.sh' for reading: $!", undef);
        while (<F>) {
            $numJobs++   if (m/^\s+job=/);
        }
        close(F);

        print STDERR "-- Configured $numJobs maln precompute jobs.\n";
    }

    if (-e "$path/maln.sh") {
        my $numJobs = 0;
        open(F, "< $path/maln.sh") or caFailure("can't open '$path/maln.sh' for reading: $!", undef);
        while (<F>) {
            $numJobs++  if (m/^\s+qry=/);
        }
        close(F);

        print STDERR "-- Configured $numJobs maln overlap jobs.\n";
    }

    stashFile("$path/precompute.sh");
    stashFile("$path/mhap.sh");

  finishStage:
    emitStage($asm, "$tag-malnConfigure");
    buildHTML($asm, $tag);

  allDone:
    stopAfter("overlapConfigure");
}


sub malnPrecomputeCheck ($$$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $typ     = shift @_;
    my $attempt = getGlobal("canuIteration");

    my $base;                #  e.g., $base/1-overlapper/mhap.sh
    my $path;                #  e.g., $path/mhap.sh

    $base = "correction"  if ($tag eq "cor");
    $base = "trimming"    if ($tag eq "obt");
    $base = "unitigging"  if ($tag eq "utg");

    $path = "$base/1-overlapper";

    goto allDone   if (skipStage($asm, "$tag-malnPrecomputeCheck", $attempt) == 1);
    goto allDone   if (fileExists("$path/precompute.files"));
    goto allDone   if (-e "$base/$asm.ovlStore");
    goto allDone   if (fileExists("$base/$asm.ovlStore.tar"));

    fetchFile("$path/precompute.sh");

    #  Figure out if all the tasks finished correctly.

    my $currentJobID   = 1;
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    open(F, "< $path/precompute.sh") or caFailure("can't open '$path/precompute.sh' for reading: $!", undef);
    while (<F>) {
        if (m/^\s+job=\"(\d+)\"$/) {
            if (fileExists("$path/blocks/$1.fasta")) {
                push @successJobs, "1-overlapper/blocks/$1.fasta\n";
            } else {
                $failureMessage .= "--   job 1-overlapper/blocks/$1.fasta FAILED.\n";
                push @failedJobs, $currentJobID;
            }

            $currentJobID++;
        }
    }
    close(F);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If not the first attempt, report the jobs that failed, and that we're recomputing.

        if ($attempt > 1) {
            print STDERR "--\n";
            print STDERR "-- ", scalar(@failedJobs), " maln precompute jobs failed:\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  If too many attempts, give up.

        if ($attempt > getGlobal("canuIterationMax")) {
            caExit("failed to precompute maln indices.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
        }

        #  Otherwise, run some jobs.

        emitStage($asm, "$tag-malnPrecomputeCheck", $attempt);
        buildHTML($asm, $tag);

        submitOrRunParallelJob($asm, "${tag}maln", $path, "precompute", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- All ", scalar(@successJobs), " maln precompute jobs finished successfully.\n";

    open(L, "> $path/precompute.files") or caExit("failed to open '$path/precompute.files'", undef);
    print L @successJobs;
    close(L);

    stashFile("$path/precompute.files");

    emitStage($asm, "$tag-malnPrecomputeCheck");
    buildHTML($asm, $tag);

  allDone:
}



sub malnCheck ($$$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $typ     = shift @_;
    my $attempt = getGlobal("canuIteration");

    my $base;                #  e.g., $base/1-overlapper/mhap.sh
    my $path;                #  e.g., $path/mhap.sh

    $base = "correction"  if ($tag eq "cor");
    $base = "trimming"    if ($tag eq "obt");
    $base = "unitigging"  if ($tag eq "utg");

    $path = "$base/1-overlapper";

    goto allDone   if (skipStage($asm, "$tag-malnCheck", $attempt) == 1);
    goto allDone   if (fileExists("$path/maln.files"));
    goto allDone   if (-e "$base/$asm.ovlStore");
    goto allDone   if (fileExists("$base/$asm.ovlStore.tar"));

    fetchFile("$path/maln.sh");

    #  Figure out if all the tasks finished correctly.

    my $currentJobID   = 1;
    my @malnJobs;
    my @successJobs;
    my @miscJobs;
    my @failedJobs;
    my $failureMessage = "";

    open(F, "< $path/maln.sh") or caExit("failed to open '$path/maln.sh'", undef);
    while (<F>) {
        if (m/^\s+qry=\"(\d+)\"$/) {
            if      (fileExists("$path/results/$1.ovb.gz")) {
                push @malnJobs,    "1-overlapper/results/$1.maln\n";
                push @successJobs, "1-overlapper/results/$1.ovb.gz\n";
                push @miscJobs,    "1-overlapper/results/$1.stats\n";
                push @miscJobs,    "1-overlapper/results/$1.counts\n";

            } elsif (fileExists("$path/results/$1.ovb")) {
                push @malnJobs,    "1-overlapper/results/$1.maln\n";
                push @successJobs, "1-overlapper/results/$1.ovb\n";
                push @miscJobs,    "1-overlapper/results/$1.stats\n";
                push @miscJobs,    "1-overlapper/results/$1.counts\n";

            } elsif (fileExists("$path/results/$1.ovb.bz2")) {
                push @malnJobs,    "1-overlapper/results/$1.maln\n";
                push @successJobs, "1-overlapper/results/$1.ovb.bz2\n";
                push @miscJobs,    "1-overlapper/results/$1.stats\n";
                push @miscJobs,    "1-overlapper/results/$1.counts\n";

            } elsif (fileExists("$path/results/$1.ovb.xz")) {
                push @malnJobs,    "1-overlapper/results/$1.maln\n";
                push @successJobs, "1-overlapper/results/$1.ovb.xz\n";
                push @miscJobs,    "1-overlapper/results/$1.stats\n";
                push @miscJobs,    "1-overlapper/results/$1.counts\n";

            } else {
                $failureMessage .= "--   job 1-overlapper/results/$1.ovb FAILED.\n";
                push @failedJobs, $currentJobID;
            }

            $currentJobID++;
        }
    }
    close(F);

    #  Also find the queries symlinks so we can remove those.  And the query directories, because
    #  the last directory can be empty, and so we'd never see it at all if only finding files.

    open(F, "find $path/queries -print |");
    while (<F>) {
        push @malnJobs, $_;
    }
    close(F);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If not the first attempt, report the jobs that failed, and that we're recomputing.

        if ($attempt > 1) {
            print STDERR "--\n";
            print STDERR "-- ", scalar(@failedJobs), " maln jobs failed:\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  If too many attempts, give up.

        if ($attempt > getGlobal("canuIterationMax")) {
            caExit("failed to compute maln overlaps.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
        }

        #  Otherwise, run some jobs.

        emitStage($asm, "$tag-malnCheck", $attempt);
        buildHTML($asm, $tag);
        submitOrRunParallelJob($asm, "${tag}maln", $path, "maln", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " maln overlap output files.\n";

    open(L, "> $path/maln.files") or caExit("failed to open '$path/maln.files'", undef);
    print L @malnJobs;
    close(L);

    open(L, "> $path/ovljob.files") or caExit("failed to open '$path/ovljob.files'", undef);
    print L @successJobs;
    close(L);

    open(L, "> $path/ovljob.more.files") or caExit("failed to open '$path/ovljob.more.files'", undef);
    print L @miscJobs;
    close(L);

    stashFile("$path/maln.files");
    stashFile("$path/ovljob.files");
    stashFile("$path/ovljob.more.files");

    emitStage($asm, "$tag-malnCheck");
    buildHTML($asm, $tag);

  allDone:
    stopAfter("overlap");
}

