include Makefile.macos

TARGET  := skm
INCDIR := inc
BINDIR := bin
OBJDIR := tmp

CXXFLAGS += -O3 $(addprefix -I, $(INCDIR))

SOURCES := $(wildcard \
	src/*.cpp \
	src/communicators/*.cpp \
	src/problem/*.cpp \
	src/utils/*.cpp \
)

TARGET := $(addprefix $(BINDIR)/, $(TARGET))
OBJECTS := $(addprefix $(OBJDIR)/, $(patsubst %.cpp, %.o, $(SOURCES)))
OBJECTS_OMP := $(addprefix $(OBJDIR)/, $(patsubst %.cpp, %_omp.o, $(SOURCES)))

all: $(TARGET) $(TARGET)_omp


$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(OBJDIR)/%_omp.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) -c $(CXXFLAGS) -DCMC_OMP $< -o $@

$(TARGET): $(OBJECTS)
	@mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) -o $(TARGET)

$(TARGET)_omp: $(OBJECTS_OMP)
	@mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS_OMP) $(OMPFLAGS) -o $(TARGET)_omp

clean:
	rm -rf $(OBJDIR) $(BINDIR)
