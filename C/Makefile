CFLAGS = -O3 -I.

RM = rm -f

CC = gcc

OBJS = egsa87.o egsa87_test.o
SRCS = egsa87.c egsa87_test.c

egsa87_test: ${OBJS}
	${CC} ${LDFLAGS} -o egsa87_test ${OBJS}

clean:
	${RM} ${OBJS} egsa87_test
