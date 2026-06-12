PROJECT = shelly-dimmer-stm32
BUILD_DIR = bin

CFILES = main.c dimcomp.c

BINARY = usart

# dimcomp.c pulls in sinf/cosf/fmodf from the math library.
LDLIBS += -lm

# TODO - you will need to edit these two lines!
DEVICE=stm32f031k6
# DEVICE=stm32f051k8

# You shouldn't have to edit anything below here.
OPENCM3_DIR=../libopencm3

include $(OPENCM3_DIR)/mk/genlink-config.mk
include ../rules.mk
include $(OPENCM3_DIR)/mk/genlink-rules.mk
