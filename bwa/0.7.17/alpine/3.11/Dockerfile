FROM alpine:3.11 as builder
RUN apk add --no-cache build-base zlib-dev tar curl
WORKDIR builddir
RUN set -ex \
    && curl -L https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 | tar xjf - --strip-components=1 \
    && make -j

# picard is included to enable piping from SamToFastq to handle uBAMs (also
# MergeBamAlignments see GATK best practices).

FROM alpine:3.11
RUN apk add --no-cache openjdk8-jre bash
ADD https://github.com/broadinstitute/picard/releases/download/2.21.6/picard.jar /picard.jar
COPY --from=builder /builddir/bwa /usr/local/bin
