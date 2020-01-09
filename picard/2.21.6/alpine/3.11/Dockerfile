FROM alpine:3.11
RUN apk add --no-cache openjdk8-jre
ADD https://github.com/broadinstitute/picard/releases/download/2.21.6/picard.jar /picard.jar
ADD picard.sh /usr/local/bin/picard
