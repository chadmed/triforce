# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2024 James Calligeros <jcalligeros99@gmail.com>

PREFIX ?= /usr/local
LIBDIR ?= $(PREFIX)/lib64

default:
	cargo build --release

performance-test:
	cargo build --release --example performance_test
	hyperfine -w 3 -- 'target/release/examples/performance_test 600'

install:
	install -dDm0755 $(DESTDIR)/$(LIBDIR)/lv2/triforce.lv2/
	install -pm0755 target/release/libtriforce.so $(DESTDIR)/$(LIBDIR)/lv2/triforce.lv2/triforce.so
	install -pm0644 triforce.ttl $(DESTDIR)/$(LIBDIR)/lv2/triforce.lv2/triforce.ttl
	install -pm0644 manifest.ttl $(DESTDIR)/$(LIBDIR)/lv2/triforce.lv2/manifest.ttl

uninstall:
	rm -rf $(DESTDIR)/$(LIBDIR)/lv2/triforce.lv2/
