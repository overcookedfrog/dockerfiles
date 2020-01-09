FROM alpine:3.11 as builder
RUN apk add --no-cache build-base curl zlib-dev bzip2-dev xz-dev bash python2
WORKDIR /app
RUN curl -L https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz | tar xzf - --strip-components=1 \
    && make -j

FROM alpine:3.11
RUN apk add --no-cache libstdc++ zlib libbz2 xz-libs
COPY --from=builder /app/bin/* /usr/local/bin/
