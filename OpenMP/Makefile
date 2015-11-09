# ===
# Rules

run: driver.c RngStream.c
	R CMD SHLIB $^

bar: bar.c
	R CMD SHLIB $^

.PHONY: clean

clean:
	rm -f *.o *.so bar.o* *~