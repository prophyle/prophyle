MAKE_PID := $(shell echo $$PPID)
JOB_FLAG := $(filter -j%, $(subst -j ,-j,$(shell ps T | grep "^\s*$(MAKE_PID).*$(MAKE)")))
JOBS     := $(subst -j,,$(JOB_FLAG))
ifeq ($(JOBS),)
	JOB_FLAG := $(filter --jobs%, $(subst --jobs ,--jobs,$(shell ps T | grep "^\s*$(MAKE_PID).*$(MAKE)")))
	JOBS     := $(subst --jobs,,$(JOB_FLAG))
endif
ifeq ($(JOBS),)
	JOBS := 1
endif
