JFLAGS = -g
JC = javac
RM = rm
.SUFFIXES: .java .class
.java.class:
	$(JC) $(JFLAGS) -cp Jama.jar -g $*.java

CLASSES = \
		  PETS.java \

default: classes
classes: $(CLASSES:.java=.class)

clean:
	$(RM) *.class
