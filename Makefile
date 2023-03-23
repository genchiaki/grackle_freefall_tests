EXEC   = star

SRCS   = onezone.c 
 


OBJS   = $(SRCS:.c=.o)
#INCL   = allvars.h nr.h


.KEEP_STATE:

CFLAGS =  $(OPT1) -O3 -DCONFIG_BFLOAT_8 -I$(HOME)/usr/local/include

LIBS   =  -lm -lgrackle -lhdf5 -L$(HOME)/usr/local/lib

CC     =  gcc

%.o : %.C
	$(CC) $(CFLAGS) -c $*.C


$(EXEC): $(OBJS) 
	gcc $(CFLAGS) $(OBJS)  -o $(EXEC)  $(LIBS)  $(OPT)

$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS) $(EXEC)
