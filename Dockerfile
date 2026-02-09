# syntax=docker/dockerfile:1.6

# ========= builder =========
FROM 192.168.3.38:5000/algo/adacus_dev:1.11.0 AS builder

ENV DEBIAN_FRONTEND=noninteractive \
    CARGO_HOME=/root/.cargo \
    RUSTUP_HOME=/root/.rustup

WORKDIR /app

# 1) 只复制构建清单，最大化缓存命中
COPY Cargo.toml Cargo.lock ./
# 如果是 workspace，把各成员 Cargo.toml 也先拷进来（否则依赖图会变，缓存不稳）
# COPY crates/*/Cargo.toml crates/*/

# 2) 预热依赖（使用 BuildKit cache）
RUN --mount=type=cache,target=/root/.cargo/registry \
    --mount=type=cache,target=/root/.cargo/git \
    --mount=type=cache,target=/app/target \
    bash -lc 'set -eux; \
    mkdir -p src && echo "fn main(){}" > src/main.rs; \
    . "$HOME/.cargo/env"; \
    cargo build --release; \
    rm -rf src'

# 3) 再复制真实源码（这一步变更最频繁，放后面）
COPY . .

# 4) 真正构建与安装（继续用 cache）
RUN --mount=type=cache,target=/root/.cargo/registry \
    --mount=type=cache,target=/root/.cargo/git \
    --mount=type=cache,target=/app/target \
    bash -lc 'set -eux; \
    . "$HOME/.cargo/env"; \
    make build-and-install; \
    strip -s /usr/bin/primer_demux || true'


# ========= runtime =========
# 方案A：更小的 Debian slim（推荐先试这个）
FROM debian:bookworm-slim AS runtime

ENV RUST_LOG=info \
    RUST_BACKTRACE=1

# 如果你的程序依赖 openssl/zlib 等动态库，在这里补齐
#（根据 ldd 输出调整）
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# 非 root 运行（建议）
RUN useradd -m -u 10001 appuser
WORKDIR /home/appuser

COPY --from=builder /usr/bin/primer_demux /usr/bin/primer_demux

USER appuser
ENTRYPOINT ["/usr/bin/primer_demux"]