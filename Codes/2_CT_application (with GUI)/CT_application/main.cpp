#include "CT_application.h"
#include "stdafx.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    CT_application w;
    w.resize(1000, 800);
    w.show();


	



    return a.exec();
}
