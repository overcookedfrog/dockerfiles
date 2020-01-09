FROM alpine:3.11 as builder
WORKDIR /app
ADD GThreads.patch .
RUN apk add --no-cache build-base zlib-dev curl
RUN curl -L http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.0.6.tar.gz | tar xzf - --strip-components=1 \
    && patch -p 1 < GThreads.patch \
    && make release

FROM alpine:3.11
RUN apk add --no-cache libstdc++ bash
COPY --from=builder /app/stringtie /usr/bin/stringtie

