# Bioinformatics-and-Server-Extras


Extra scripts and notes to get people going with poolseq, as well as various handy scripts for server maintenance

## `sendipcheck.sh:`

> A shell script you can set as a cronjob  that will check if your IP has changed at a given interval, then send that to a Slack channel you can configure with a Slack webhook (see: [Sending messages using Incoming Webhooks](https://api.slack.com/messaging/webhooks)). 

## `workflow_single_library_assemblies.sh:`

> THis wrapper is designed to start from a given `<project>` directory, this script will 1) take SE or PE reads from `./<project>/raw/` (with file names containing `R1` for F reads and `R2` for reverse), 2) trim them using [Trim Galore!](http://quast.sourceforge.net/quast), 3) run a *de novo* assembly for each library (SE or PE) using [SPAdes](http://cab.spbu.ru/software/spades/) with both default and user selected kmer lengths, and 4) produce summary statistics that can be viewed in an interactive HMTL output produced by [QUAST](http://bioinf.spbau.ru/quast). This has been batch tested on over 80 microbial genomes in a run, as well as a dozen *de novo* assemblies of marine inverebrates in a run. 

*note: next upload will have an alternative workflow for multilibray (in application where multiple libraries contribute to the desired de novo assembly, as well as further steps for mapping, and variant calling*
