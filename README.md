# Bioinformatics-and-Server-Extras


Extra scripts and notes to get people going with poolseq, as well as various handy scripts for server maintenance

## `sendipcheck.sh:`

> A shell script you can set as a cronjob  that will check if your IP has changed at a given interval, then send that to a Slack channel you can configure with a Slack webhook (see: [Sending messages using Incoming Webhooks](https://api.slack.com/messaging/webhooks)). 