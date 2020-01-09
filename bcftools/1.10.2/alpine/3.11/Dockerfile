FROM alpine:3.11 as builder
RUN apk add --no-cache build-base zlib-dev bzip2-dev xz-dev curl-dev curl
WORKDIR /app
RUN curl -L https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 | tar xj --strip-components=1 \
    && ./configure \
    && make -j \
    && make install

FROM alpine:3.11
RUN apk add --no-cache zlib libbz2 xz-libs libcurl bash
COPY --from=builder /usr/local/bin/bcftools /usr/local/lib/bcftools
