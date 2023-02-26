CC := gcc
CFLAGS = -Wall
SRCS := $(wildcard *.c)
OBJS := $(SRCS:%.c=%.o)
TARGET := simpleH264Analyzer

all: $(TARGET)

$(TARGET):$(OBJS)
	$(CC) $(CFLAGS) $^ -o $@

.c.o:
	$(CC)  $(CFLAGS) -c -MMD -o $@ $<

.PHONY: clean

clean:
	rm -fr *.o *.d $(TARGET) trace.txt
