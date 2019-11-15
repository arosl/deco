CFLAGS := -Wall -Werror -O3
LDFLAGS := -lm

TARGETS := deco

all: $(TARGETS)

%.o: %.c
	gcc $(CFLAGS) -c $<

deco: deco.o
	gcc $(LDFLAGS) -o $@ $^

clean:
	rm -f *~ *.o

proper: clean
	rm -f $(TARGETS)
