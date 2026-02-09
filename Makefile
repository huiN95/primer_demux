build:
	cargo build --release

install:
	cp target/release/primer_demux /usr/bin/

build-and-install:
	cargo build --release
	cp target/release/primer_demux /usr/bin/
