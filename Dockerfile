
# ========= builder =========
FROM 192.168.3.38:5000/algo/adacus_dev:1.11.0 AS builder

WORKDIR /app
COPY . .

RUN bash -lc 'set -eux; \
    . "$HOME/.cargo/env"; \
    make build-and-install; \
    strip -s /usr/bin/primer_demux || true'

# ========= runtime =========
FROM debian:bookworm-slim AS runtime

ENV RUST_LOG=info \
    RUST_BACKTRACE=1

# Install runtime dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user and work directory
RUN useradd -m -u 10001 -s /usr/sbin/nologin appuser \
    && mkdir -p /work \
    && chown -R appuser:appuser /work

WORKDIR /work

COPY --from=builder /usr/bin/primer_demux /usr/bin/primer_demux

USER appuser
ENTRYPOINT ["/usr/bin/primer_demux"]
