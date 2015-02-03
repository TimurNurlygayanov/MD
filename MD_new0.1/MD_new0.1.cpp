// MD_new0.1.cpp: главный файл проекта.

#include "prog.cpp"
#include "stdafx.h"

using namespace System;

int main(array<System::String ^> ^args)
{
	K = 12;
	K2 = 46;
	load_seed("save_file.txt");
	init("program.txt");

    return 0;
}
