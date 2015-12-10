# ===
# Rules

bar: bar.c RngStream.c
	R CMD SHLIB $^

.PHONY: clean

clean:
	rm -f *.o *.so bar.o* *~