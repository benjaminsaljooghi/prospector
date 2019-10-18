#
# Makefile:  /home/benjamin/foo/Makefile.blast_demo.app
#
# This file was originally generated by shell script "new_project.sh" (r565247)
# Fri Aug 16 11:32:33 AEST 2019
#

###  BASIC PROJECT SETTINGS
APP = blast_demo.out
SRC = blast_demo
OBJ_DIR = ../obj
OBJ = $(OBJ_DIR)/dlinked $(OBJ_DIR)/prospector $(OBJ_DIR)/util
MY_LIBS = -L/usr/local/cuda/lib64 -lcudart


LIB_ = $(BLAST_INPUT_LIBS) $(BLAST_LIBS) $(OBJMGR_LIBS)
LIB = $(LIB_:%=%$(STATIC))
LIBS = $(BLAST_THIRD_PARTY_LIBS) $(CMPRS_LIBS) $(NETWORK_LIBS) $(DL_LIBS) $(ORIG_LIBS) $(MY_LIBS)

# These settings are necessary for optimized WorkShop builds, due to
# BLAST's own use of them.
CXXFLAGS = $(FAST_CXXFLAGS)

### LOCAL_LDFLAGS automatically added
LDFLAGS = $(LOCAL_LDFLAGS) $(FAST_LDFLAGS)

REQUIRES = objects -Cygwin

###  EXAMPLES OF OTHER SETTINGS THAT MIGHT BE OF INTEREST
# PRE_LIBS = $(NCBI_C_LIBPATH) .....
# CFLAGS   = $(FAST_CFLAGS)
# CXXFLAGS = $(FAST_CXXFLAGS)
# LDFLAGS  = $(FAST_LDFLAGS)