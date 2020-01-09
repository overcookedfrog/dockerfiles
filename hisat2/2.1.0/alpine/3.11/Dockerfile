FROM alpine:3.11 as builder
RUN apk add --update build-base curl zlib-dev bzip2-dev ncurses-dev xz-dev curl-dev
WORKDIR /samtools-1.9
RUN curl -L https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 | tar xj --strip-components=1 \
    && ./configure --prefix=/app/samtools-1.9 \
    && make -j
WORKDIR /hisat2
RUN curl -LO ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-source.zip && \
    unzip hisat2-2.1.0-source.zip && \
    cd hisat2-2.1.0 && \
    make -j

FROM alpine:3.11
RUN apk add --no-cache procps libstdc++ zlib libcurl libbz2 xz-libs ncurses bash perl python3
COPY --from=builder /samtools-1.9/samtools /usr/local/bin/samtools
COPY --from=builder /hisat2/hisat2-2.1.0/extract_exons.py                       /usr/local/bin/extract_exons.py
COPY --from=builder /hisat2/hisat2-2.1.0/extract_splice_sites.py                /usr/local/bin/extract_splice_sites.py
COPY --from=builder /hisat2/hisat2-2.1.0/hisat2                                 /usr/local/bin/hisat2
COPY --from=builder /hisat2/hisat2-2.1.0/hisat2-align-l                         /usr/local/bin/hisat2-align-l
COPY --from=builder /hisat2/hisat2-2.1.0/hisat2-align-s                         /usr/local/bin/hisat2-align-s
COPY --from=builder /hisat2/hisat2-2.1.0/hisat2-build                           /usr/local/bin/hisat2-build
COPY --from=builder /hisat2/hisat2-2.1.0/hisat2-build-l                         /usr/local/bin/hisat2-build-l
COPY --from=builder /hisat2/hisat2-2.1.0/hisat2-build-s                         /usr/local/bin/hisat2-build-s
COPY --from=builder /hisat2/hisat2-2.1.0/hisat2_extract_exons.py                /usr/local/bin/hisat2_extract_exons.py
COPY --from=builder /hisat2/hisat2-2.1.0/hisat2_extract_snps_haplotypes_UCSC.py /usr/local/bin/hisat2_extract_snps_haplotypes_UCSC.py
COPY --from=builder /hisat2/hisat2-2.1.0/hisat2_extract_snps_haplotypes_VCF.py  /usr/local/bin/hisat2_extract_snps_haplotypes_VCF.py
COPY --from=builder /hisat2/hisat2-2.1.0/hisat2_extract_splice_sites.py         /usr/local/bin/hisat2_extract_splice_sites.py
COPY --from=builder /hisat2/hisat2-2.1.0/hisat2-inspect                         /usr/local/bin/hisat2-inspect
COPY --from=builder /hisat2/hisat2-2.1.0/hisat2-inspect-l                       /usr/local/bin/hisat2-inspect-l
COPY --from=builder /hisat2/hisat2-2.1.0/hisat2-inspect-s                       /usr/local/bin/hisat2-inspect-s
COPY --from=builder /hisat2/hisat2-2.1.0/hisat2_simulate_reads.py               /usr/local/bin/hisat2_simulate_reads.py
COPY --from=builder /hisat2/hisat2-2.1.0/hisatgenotype_build_genome.py          /usr/local/bin/hisatgenotype_build_genome.py
COPY --from=builder /hisat2/hisat2-2.1.0/hisatgenotype_extract_reads.py         /usr/local/bin/hisatgenotype_extract_reads.py
COPY --from=builder /hisat2/hisat2-2.1.0/hisatgenotype_extract_vars.py          /usr/local/bin/hisatgenotype_extract_vars.py
COPY --from=builder /hisat2/hisat2-2.1.0/hisatgenotype_hla_cyp.py               /usr/local/bin/hisatgenotype_hla_cyp.py
COPY --from=builder /hisat2/hisat2-2.1.0/hisatgenotype_locus.py                 /usr/local/bin/hisatgenotype_locus.py
ADD hisat2_ubam.py /usr/local/bin/hisat2-ubam
ADD hisat2_ubams.py /usr/local/bin/hisat2-ubams
