
TEMPLATE	= app
LANGUAGE	= C++

CONFIG	+= warn_on debug

INCLUDEPATH	+= $(LITTLEDOG)/include $(LITTLEDOG)/rt/robots/littledog/api/lib ../reference_ui

SOURCES	+= main.cpp

unix {
	UI_DIR = .ui
	MOC_DIR = .moc
	OBJECTS_DIR = .obj

	exists($(LITTLEDOG)/lib/linux/x86-64_gcc34/liblittledog.a) {
		LIBS += $(LITTLEDOG)/lib/linux/x86-64_gcc34/liblittledog.a
	}

	exists($(LITTLEDOG)/lib/linux/i586_gcc34/liblittledog.a) {
		LIBS += $(LITTLEDOG)/lib/linux/i586_gcc34/liblittledog.a
	}
}

