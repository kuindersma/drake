
TEMPLATE	= app
LANGUAGE	= C++

CONFIG	+= qt warn_on debug

INCLUDEPATH	+= $(LITTLEDOG)/include $(LITTLEDOG)/rt/robots/littledog/api/lib ../reference_ui

HEADERS	+= ../reference_ui/LittleDogUI.h \
	MyLittleDog.h \
	MyLittleDogUI.h

SOURCES	+= main.cpp \
	../reference_ui/LittleDogUI.cpp \
	MyLittleDog.cpp \
	MyLittleDogUI.cpp \
	KinematicWalker.cpp \
	SwingStanceGenerator.cpp \
	LinearSpline.cpp \

FORMS = ../reference_ui/LittleDogUI_uic.ui \
	../reference_ui/DialogCalibratePose_uic.ui

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
